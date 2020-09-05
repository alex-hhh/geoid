#lang info
(define collection "geoid")
(define deps '("base" "math-lib" "rackunit-lib" "typed-racket-lib"))
(define build-deps '("racket-doc" "scribble-lib" "al2-test-runner"))
(define scribblings '(("scribblings/geoid.scrbl" ())))
(define pkg-desc "work efficiently with geographic data")
(define version "0.0")
(define pkg-authors '(AlexHarsanyi@gmail.com))
