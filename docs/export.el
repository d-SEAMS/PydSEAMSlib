;; Batch export org-mode files to RST for Sphinx.
;; Usage (cwd = docs/): emacs --batch --load export.el
;; Org under orgmode/ is the source. RST under source/ is generated and committed.
;; Pages CI builds the committed RST; it does not run this file.
;; Exclude reference/api.org -- the autodoc landing page is source/api.md (Myst).
(require 'package)
(add-to-list 'package-archives '("melpa" . "https://melpa.org/packages/") t)
(package-initialize)

;; ox-rst 2025-04 calls the Org 9.7 element API (org-element-parent,
;; org-element-type-p, ...). Emacs 29 bundles Org 9.6, so take Org from
;; GNU ELPA before anything loads the bundled copy.
(unless (package-installed-p 'org '(9 7))
  (package-refresh-contents)
  (package-install (cadr (assq 'org package-archive-contents))))
(unless (package-installed-p 'ox-rst)
  (unless package-archive-contents
    (package-refresh-contents))
  (package-install 'ox-rst))

(require 'ox-rst)
(require 'ox-publish)
(unless (version<= "9.7" org-version)
  (error "ox-rst needs Org 9.7 or newer; loaded Org %s" org-version))

(setq org-export-with-section-numbers nil)
(setq org-export-with-toc nil)
(setq org-export-with-author nil)
(setq org-export-with-timestamps nil)
(setq org-rst-headline-underline ?-)

(setq org-publish-project-alist
      '(("pydseams-rst"
         :base-directory "./orgmode/"
         :base-extension "org"
         :publishing-directory "./source/"
         :publishing-function org-rst-publish-to-rst
         :recursive t
         :exclude "^reference/api\\.org$"
         :headline-levels 4
         :with-toc nil
         :section-numbers nil
         :with-author nil)))

(org-publish "pydseams-rst" t)

;; Verify no Org syntax leaked into the generated reStructuredText.
(let ((leaks '()))
  (dolist (file (directory-files-recursively "./source" "\\.rst$"))
    (with-temp-buffer
      (insert-file-contents file)
      (goto-char (point-min))
      (while (re-search-forward "^#\\+\\(begin\\|end\\)_\\|^#\\+[A-Z_]+:\\|:ref:\\\\`\\|\\[\\[file:" nil t)
        (push (format "%s:%d: %s"
                      file
                      (line-number-at-pos)
                      (buffer-substring (line-beginning-position) (line-end-position)))
              leaks))))
  (when leaks
    (error "Org syntax leaked into generated RST:\n%s"
           (mapconcat #'identity (nreverse leaks) "\n"))))
