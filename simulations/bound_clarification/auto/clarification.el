(TeX-add-style-hook "clarification"
 (lambda ()
    (LaTeX-add-environments
     "theorem")
    (TeX-add-symbols
     "E")
    (TeX-run-style-hooks
     "graphicx"
     "amsthm"
     "amsmath"
     "fullpage"
     "latex2e"
     "amsart10"
     "amsart")))

