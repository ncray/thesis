(TeX-add-style-hook "friedman-test"
 (lambda ()
    (LaTeX-add-labels
     "C:friedman-test"
     "item:friedman-score"
     "D:kernel"
     "S:SVM-section"
     "eq:l1-svm-primal"
     "eq:l1-svm-dual"
     "eq:hinge-svm"
     "T:representer"
     "eq:representer"
     "friedman_equiv"
     "S:MMD"
     "fig:c_param"
     "fig:null_dist"
     "fig:power_normal"
     "twitter_data"
     "fig:power_string"
     "fig:birds"
     "fig:power_birds")
    (TeX-run-style-hooks
     "friedman-test/img/c_param"
     "friedman-test/img/null_dist"
     "friedman-test/img/power_normal"
     "friedman-test/img/power_string"
     "friedman-test/img/power_birds")))

