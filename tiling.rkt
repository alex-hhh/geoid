#lang typed/racket/base

;; SPDX-License-Identifier: LGPL-3.0-or-later
;;
;; tiling.rkt -- functions to determine the list of geoids which cover regions
;; on the Earth surface.
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

(require "private/tiling.rkt")

(provide
 ;; this is exported so contracts can be written by code using this library,
 ;; there are no user-invocable operations on the class
 region%

 make-spherical-cap
 make-closed-polyline
 make-open-polyline

 subtract-regions
 intersect-regions
 join-regions

 guess-winding-order

 geoid-tiling-for-region
 refine-geoid-tiling-for-region
 coarse-geoid-tiling-for-region)
