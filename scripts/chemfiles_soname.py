"""Print the vendored chemfiles library and the name the loader opens.

seams-core builds patched chemfiles into chemfiles-prefix. The file
the dynamic loader asks for is the ELF SONAME or the Mach-O install id,
which is not the same string on Linux and macOS. stdout is two lines,
the real file and that loader name, or empty when the prefix has no
library.
"""

from __future__ import annotations

import struct
import sys
from pathlib import Path


def _elf_soname(data: bytes) -> str | None:
    if data[:4] != b"\x7fELF":
        return None
    elf_class = data[4]
    endian = "<" if data[5] == 1 else ">"
    if elf_class == 2:
        e_phoff = struct.unpack_from(endian + "Q", data, 32)[0]
        e_phentsize = struct.unpack_from(endian + "H", data, 54)[0]
        e_phnum = struct.unpack_from(endian + "H", data, 56)[0]
        # p_type, p_flags, p_offset, p_vaddr, p_paddr, p_filesz, p_memsz, p_align
        ph_fmt = endian + "IIQQQQQQ"
        offset_i, vaddr_i, filesz_i = 2, 3, 5
        ent = 16
        tag = endian + "Q"
    elif elf_class == 1:
        e_phoff = struct.unpack_from(endian + "I", data, 28)[0]
        e_phentsize = struct.unpack_from(endian + "H", data, 42)[0]
        e_phnum = struct.unpack_from(endian + "H", data, 44)[0]
        ph_fmt = endian + "IIIIIIII"
        offset_i, vaddr_i, filesz_i = 1, 2, 4
        ent = 8
        tag = endian + "I"
    else:
        return None
    loads: list[tuple[int, int, int]] = []
    dyn_off = dyn_size = None
    for index in range(e_phnum):
        at = e_phoff + index * e_phentsize
        fields = struct.unpack_from(ph_fmt, data, at)
        if fields[0] == 1:
            loads.append((fields[vaddr_i], fields[offset_i], fields[filesz_i]))
        elif fields[0] == 2:
            dyn_off = fields[offset_i]
            dyn_size = fields[filesz_i]
    if dyn_off is None or dyn_size is None:
        return None

    def file_offset(vaddr: int) -> int | None:
        for base, start, size in loads:
            if base <= vaddr < base + size:
                return start + (vaddr - base)
        return None

    strtab = soname_at = None
    pos = dyn_off
    end = dyn_off + dyn_size
    while pos + ent <= end:
        kind = struct.unpack_from(tag, data, pos)[0]
        value = struct.unpack_from(tag, data, pos + ent // 2)[0]
        if kind == 5:
            strtab = value
        elif kind == 14:
            soname_at = value
        elif kind == 0:
            break
        pos += ent
    if strtab is None or soname_at is None:
        return None
    start = file_offset(strtab)
    if start is None:
        return None
    blob = data[start + soname_at :]
    return blob.split(b"\x00", 1)[0].decode("ascii")


def _macho_id(data: bytes) -> str | None:
    magic = struct.unpack_from("<I", data, 0)[0]
    if magic == 0xBEBAFECA:
        count = struct.unpack_from(">I", data, 4)[0]
        if count < 1:
            return None
        offset = struct.unpack_from(">I", data, 16)[0]
        return _macho_id(data[offset:])
    if magic == 0xFEEDFACF:
        header = 32
        endian = "<"
    elif magic == 0xFEEDFACE:
        header = 28
        endian = "<"
    elif magic == 0xCFFAEDFE:
        header = 32
        endian = ">"
    elif magic == 0xCEFAEDFE:
        header = 28
        endian = ">"
    else:
        return None
    ncmds = struct.unpack_from(endian + "I", data, 16)[0]
    pos = header
    for _ in range(ncmds):
        cmd, cmdsize = struct.unpack_from(endian + "II", data, pos)
        if cmd == 0xD:  # LC_ID_DYLIB
            name_off = struct.unpack_from(endian + "I", data, pos + 8)[0]
            blob = data[pos + name_off :]
            path = blob.split(b"\x00", 1)[0].decode("utf-8", "replace")
            return Path(path).name
        pos += cmdsize
    return None


def loader_name(path: Path) -> str | None:
    data = path.read_bytes()
    return _elf_soname(data) or _macho_id(data)


def real_library(libdir: Path) -> Path | None:
    found: list[Path] = []
    if not libdir.is_dir():
        return None
    for path in libdir.iterdir():
        if path.is_symlink() or not path.is_file():
            continue
        if "chemfiles" not in path.name or path.name.endswith(".a"):
            continue
        found.append(path)
    if not found:
        return None
    found.sort(key=lambda item: item.stat().st_size, reverse=True)
    return found[0]


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} LIBDIR", file=sys.stderr)
        return 2
    library = real_library(Path(sys.argv[1]))
    if library is None:
        return 1
    name = loader_name(library)
    if not name:
        return 1
    sys.stdout.write(f"{library}\n{name}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
