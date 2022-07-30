#lang typed/racket/base
;; This file is part of geoid -- work efficiently with geographic data
;; Copyright (c) 2020 Alex Harsányi <AlexHarsanyi@gmail.com>
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

(require "private/geoid.rkt")

(provide
 first-valid-geoid                      ; tested, doc
 last-valid-geoid                       ; tested, doc
 sentinel-geoid                         ; tested, doc
 valid-geoid?                           ; tested, doc
 geoid-level                            ; tested, doc
 ;; geoid-face                          ; tested, not exported
 geoid-stride                           ; tested, doc


 enclosing-geoid                        ; tested, doc
 split-geoid                            ; tested, doc
 leaf-geoid?                            ; tested, doc
 leaf-span                              ; tested, doc
 leaf-span*                             ; tested, doc
 contains-geoid?                        ; tested, doc
 leaf-corners                           ; doc

 adjacent-geoids                        ; tested, doc

 approximate-area-for-geoid             ; doc
 approximate-area-for-level             ; doc

 ;; Testing Helpers

 random-geoid                           ; tested, doc
 leaf-outline                           ; doc

 ;; Storing/Retrieving from a SQLite database

 geoid->sqlite-integer                  ; doc
 sqlite-integer->geoid                  ; doc

 ;; Conversion from latitude/longitude and back

 lat-lng->geoid                         ; tested, doc
 geoid->lat-lng                         ; tested, doc
 lat-lng-rect                           ; doc

 distance-between-geoids                ; doc
 distance-from-geoid                    ; TODO: doc
 )


;; raco setup --check-pkg-deps --pkgs geoid
;; raco test --no-run-if-absent --package geoid
