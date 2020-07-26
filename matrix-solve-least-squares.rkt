; Author: Matthew Might
;         matt.might.net

; License: CC0
;          But, if you feel like giving me a nod and a link,
;          that's certainly appreciated.

; I wrote this to get a deeper understanding of the methods used
; find least squares solutions.  The code is extensively commented
; for readability.

; TODO: Add an SVD-based approach as well.
;       Starting with a simple (if inefficient) implementation:
;
;          https://rpubs.com/aaronsc32/singular-value-decomposition-r

; TODO: Add bindings to BLAS/LAPACK for high-performance implementation.


; WARNING: I'm not suggesting using these algorithms in practice.
; I haven't paid attention to efficiency or numeral stability.
; I'm not saying anything is unstable, but I haven't even checked.


(module
 matrix-solve-least-squares
 racket/base

(provide 
 matrix-solve-least-squares/qr
 matrix-solve-least-squares/transpose)


(require math/matrix)
(require math/array)



; Given an m (rows) x n (cols) matrix A and a column vector b
; of height m, the least squares solution is a row vector x of
; length n such that the distance between A*x and b -- ∥ A*x - b ∥
; -- is minimized.

; In this case, ∥ v ∥ is the length of vector v (formally,
; the L² norm) a.k.a. the Euclidian distance.

; The L² norm of a vector <d1, d2, ..., dm> is:

;  ∥ <d1, d2, ..., dm> ∥ = √(d1² + d² + ... + dm²)

; so minimizing this norm is the same as minimizing *without* the sqrt:

;   d1² + d2² + ... + dm²

; Hence, the "least squares solution" minimizes the sum of these squares.



; Informally, the ith square in this sum is a measure of the error for
; the ith row in A and the corresponding ith entry in b.

; In fact, if you select the ith row from A (a1, a2, ..., an) and the
; corresponding row in b -- bi, then for a solution vector
; (x1, x2, ..., xn), di is the "error" in the estimation for this row:
;
;  error for point i = a1*x1 + a2*x2 + ... + an*xn - bi.



; Find the column vector x such that ∥ A*x - b ∥
; is minimized using a QR decomposition:
(define (matrix-solve-least-squares/qr A b)

  ; A is an m x n matrix:
  (define m (matrix-num-rows A))
  (define n (matrix-num-cols A))

  ;   n columns  
  ; 
  ; [           ] 
  ; [           ]
  ; [     A     ]  m rows
  ; [           ]
  ; [           ]

  ; b must be a column vector with height m.

  
  ; Construct a QR decomposition
  (define-values (Q R) (matrix-qr A))

  ; Q will be an m x m orthogonal matrix.
  ; R will be an upper triangular matrix

  ; As an orthogonal matrix, Q obeys properties that
  ; simplify the calculation of least squares.

  ; An orthogonal matrix is a square matrix where the vectors
  ; that compose the matrix are orthonormal to each other:

  ; These vectors all have length 1, and they are all orthogonal to
  ; each other.

  ; The inverse of an orthogonal matrix is its transpose:
  ;  Q⁻¹ = Qᵀ
  (define Q⁻¹ (matrix-transpose Q))

  ; The inverse of an orgothonal matrix is also orthogonal.

  ; If you multiply a vector by an orthogonal matrix, then
  ; its length is preserved:

  ; ∥ Q * v ∥ = ∥ v ∥

  
  ; We can restate the problem now as minimizing:
  ; ∥ Q * R * x - b ∥

  ; Because applying orthogonal matrices don't change the length:
  
  ;   ∥ Q * R * x - b ∥
  
  ; = ∥ Q⁻¹ * (Q * R * x - b) ∥
  
  ; = ∥ Q⁻¹ * Q * R * x  -  Q⁻¹ * b ∥
  
  ; = ∥ R * x  -  Q⁻¹ * b ∥


  ; We can simplify Q⁻¹ * b into a vector:
  (define Q⁻¹*b (matrix* Q⁻¹ b))

  ; Because R is upper triangular, finding the vector x to minimize
  ; this quantity will be an easier problem to solve now.
  
  ; Beacuse R is upper triangular, everything past row n is 0:

  ;              n columns  
  ; 
  ;     [ r11 ...                 ] 
  ;     [  0  r22 ...             ]
  ;     [  0   0  r33 ...         ]  n rows
  ;     [           ...           ]
  ; R = [  0   0   0   0  ... rnn ]
  
  ;     [  0   0   0   0  ...  0  ]
  ;     [           ...           ]  (m-n) rows
  ;     [  0   0   0   0  ...  0  ]
  

  ; Since the bottom rows are all 0s, we can define
  ; an equivalent Rʹ as an n x n (square) matrix:

  ;      [ r11 ...                 ] 
  ;      [  0  r22 ...             ]
  ; Rʹ = [  0   0  r33 ...         ]  n rows
  ;      [           ...           ]
  ;      [  0   0   0   0  ... rnn ]

  
  
  ; We can keep just the first n rows of R
  ; to get the square matrix Rʹ:
  (define Rʹ (submatrix R (:: 0 n) (::)))

  ; Ande we can similarly break apart Q⁻¹*b into 
  ; n elements and the remaining (m-n) elements


  ;           [    c1    ]
  ;           [    c2    ]
  ;           [    c3    ]  n rows
  ;           [    ...   ]
  ; Q⁻¹ * b = [    cn    ]
  
  ;           [    d1    ] 
  ;           [    ...   ] (m-n) rows
  ;           [  d(m-n)  ]
  
  

  ; The rows in the vector c portion of Q⁻¹*b line up
  ; with the rows of Rʹ
  (define c (submatrix Q⁻¹*b (:: 0 n) (::)))

  ; The remainder d will come to define the error:
  (define d (submatrix Q⁻¹*b (:: n #f) (::)))


  ; With this factoring, we can rewrite the minimization problem:

  ;   min ∥ R * x  -  Q⁻¹ * b ∥
  ;    x
  
  ; = min ∥ [ Rʹ ] * x  -  [ c ] ∥
  ;    x  ∥ [ 0  ]         [ d ] ∥

  ; = min ∥ [ Rʹ*x - c ] ∥
  ;    x  ∥ [    -d    ] ∥

  ; = ∥ d ∥   because  we can minimize  Rʹ*x - c  (to 0) with:
    
  ;           x = Rʹ⁻¹ * c

  (define Rʹ⁻¹ (matrix-inverse Rʹ))

  (define Rʹ⁻¹*c (matrix* Rʹ⁻¹ c))

  (define x-min Rʹ⁻¹*c)

  ; And thus the length of d, ∥ d ∥, is the error for the closest fit:
  ; (define err (matrix-2norm d))

  x-min)






; There is a straightforward closed-form solution for the least-squares
; problem using matrix calculus.


; Start again with the minimization problem:

;     min ∥ A*x - b ∥ 
;      x

;  =  min ∥ A*x - b ∥²
;      x

;  =  min (A*x - b)ᵀ * (A*x - b)
;      x

; At this point, we can use matrix calculus to take the derivative
; with respect to the vector x, and set it equal to 0 to solve
; for the value of x that's at the minimum:

; Let L = (A*x - b)ᵀ * (A*x - b)
;
;       =  bᵀ*b - bᵀ*A*x - xᵀ*Aᵀ*b + xᵀ*Aᵀ*A*x

; It helps to see the derivative of each term individually, which can
; be computed from the "denominator" derivative rules:

;  https://en.wikipedia.org/wiki/Matrix_calculus#Identities

;   ∂(bᵀ*b)/∂x         = 0

;   ∂(- bᵀ*A*x)/∂x     = - Aᵀ * b

;   ∂(- xᵀ*Aᵀ*b)/∂x    = - Aᵀ * b

;   ∂(xᵀ*Aᵀ*A*x)/∂x    = + 2 * Aᵀ * A * x


; Putting it all together:

;   ∂L/∂x = - 2 * Aᵀ * b  +  2 * Aᵀ * A * x


; Setting ∂L/∂x = 0 yields:

;        2 * Aᵀ * A * x   =            2 * Aᵀ * b   (1)

; =>         Aᵀ * A * x   =                Aᵀ * b   (2)

; =>                  x   =   (Aᵀ * A)⁻¹ * Aᵀ * b   (3)


; The form (2) can be used with any square matrix solver, including
; the built-in matrix-solve in Racket.

; Or, the solution can be calculated explicitly with the inverse of (Aᵀ * A).


(define (matrix-solve-least-squares/transpose A b)

  (define Aᵀ (matrix-transpose A))

  (define Aᵀ*A (matrix* Aᵀ A))

  (define Aᵀ*b (matrix* Aᵀ b))

  ; Find x such that  Aᵀ*A * x   =   Aᵀ*b:
  (define x-min (matrix-solve  Aᵀ*A  Aᵀ*b))

  x-min)




)
