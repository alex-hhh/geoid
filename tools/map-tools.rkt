#lang racket/base

;; SPDX-License-Identifier: LGPL-3.0-or-later
;;
;; map-tools.rkt -- utilities for displaying geoids and regions on a map
;;
;; This file is part of geoid -- work efficiently with geographic data
;; Copyright (c) 2022 Alex Harsányi <AlexHarsanyi@gmail.com>
;;
;; This program is free software: you can redistribute it and/or modify it
;; under the terms of the GNU Lesser General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or (at your
;; option) any later version.
;;
;; This program is distributed in the hope that it will be useful, but WITHOUT
;; ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
;; FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
;; License for more details.
;;
;; You should have received a copy of the GNU Lesser General Public License
;; along with this program.  If not, see <http://www.gnu.org/licenses/>.

(require racket/gui
         racket/class
         map-widget
         geoid/geodesy
         geoid/private/geoid
         (only-in plot/utils ->pen-color))

;; Prevent this module from being instantiated when the tests are being run,
;; since the dependencies are not installed
(module test racket/base)

(provide group-pen put-track put-cap put-geoids make-map-frame)

(define (group-pen group-name color [width 1.0] [style 'solid])
  (define c
    (if (is-a? color color%)
        color
        (match-let ([(list r g b) (->pen-color color)])
          (make-color r g b))))
  (define pen (send the-pen-list find-or-create-pen c width style))
  (lambda (the-map)
    (send the-map set-group-pen group-name pen)))

(define (put-track t [group-name #f])
  (lambda (the-map)
    (send the-map add-track t group-name)))

(define (make-circular-track lat lon radius
         #:segments [segments 24] #:closed? [closed? #t])
  (define step (/ 360.0 segments))
  (define track
    (for/list ([bearing (in-range 0.0 360.0 step)])
      (define-values (clat clon) (destination-point lat lon bearing radius))
      (vector clat clon)))
  (if closed?
      (append track (list (car track)))
      track))

(define (put-cap c [group-name #f])
  (match-define (vector x y z) (send c get-center))
  (define-values (lat lng) (unit-vector->lat-lng x y z))
  (define radius (send c get-radius))
  (put-track (make-circular-track lat lng radius) group-name))

(define (leaf-outline-as-track g)
  (for/list ([g (in-list (leaf-outline g))])
    (define-values (lat lng) (geoid->lat-lng g))
    (vector lat lng)))

(define (put-geoids gs [group-name #f])
  (define tracks (map leaf-outline-as-track gs))
  (lambda (the-map)
    (for ([t (in-list tracks)])
      (send the-map add-track t group-name))))

(define (make-map-frame items #:title [title "Map"])
  (define tl (new frame% [label title] [width 800] [height 600]))
  (define the-map (new map-widget% [parent tl]))
  (map (lambda (f) (f the-map)) items)
  (send the-map resize-to-fit)
  (send tl show #t)
  tl)
