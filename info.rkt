#lang info
(define collection "geoid")
(define license 'LGPL-3.0-or-later)
(define pkg-authors '(AlexHarsanyi@gmail.com))
(define version "0.0")
(define pkg-desc "work efficiently with geographic data")

(define deps '("base" "math-lib" "rackunit-lib" "typed-racket-lib"))
(define build-deps '("racket-doc" "scribble-lib" "al2-test-runner"))
(define scribblings '(("scribblings/geoid.scrbl" ())))

;; These files are for development only and depend on other packages...
(define compile-omit-paths '("tools/"))
(define test-omit-paths '("tools/"))
