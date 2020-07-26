#lang racket

; Author: Matthew Might
;         matt.might.net

; License: CC0
;          But, if you feel like giving me a nod and a link,
;          that's certainly appreciated.


(require math/matrix)
(require math/array)


(require "matrix-solve-least-squares.rkt")


; Test data:
(define A
 (matrix
  [[ 3 4 5 ]
   [ 6 1 2 ]
   [ 2 3 0 ]
   [ 1 1 1 ]
   [ 2 4 6 ]]))

(define b
  (->col-matrix '[1 2 3 4 5]))




(matrix-solve-least-squares/qr A b)


(matrix-solve-least-squares/transpose A b)


 
  
