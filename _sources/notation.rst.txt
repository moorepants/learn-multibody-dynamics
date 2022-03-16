========
Notation
========

This explains the notation in the book. The mathematical symbol is shown and
then an example of a variable name that we use in the code.

.. todo:: Add the notation that is in the hand drawn figures.

:math:`x`, ``x``
   Scalars are normal font.
:math:`\bar{v}`, ``v``
   Vectors are indicated with a bar.
:math:`\hat{u}`, ``uhat``
   Unit vectors are indicated with a hat.
:math:`A`, ``A``
   Reference frame :math:`A`.
:math:`\hat{a}_x,\hat{a}_y,\hat{a}_z`, ``A.x``, ``A.y``, ``A.z``
   Right handed mutually perpendicular unit vectors fixed in reference frame
   :math:`A`.
:math:`|\bar{v}|`, ``v.magnitude()``
   Magnitude of a vector; Euclidean norm (2-norm).
:math:`\bar{u} \cdot \bar{v}`, ``u.dot(v)``
   Dot product of two vectors.
:math:`\bar{u} \times \bar{v}`, ``u.cross(v)``
   Cross product of two vectors.
:math:`\bar{u} \otimes \bar{v}`, ``u.outer(v)``
   Outer product of two vectors.
:math:`\mathbf{R}`, ``R``
   Matrices are capitalized letters in bold font.
:math:`{}^A\mathbf{C}^B`, ``A_C_B``
   Direction cosine matrix relating reference frames (or rigid bodies)
   :math:`B` and :math:`A` where this relation between the right handed
   mutually perpendicular unit vectors fixed in the two reference frames that
   follow this relationship:

   .. math::

      \begin{bmatrix}
        \hat{a}_x \\
        \hat{a}_y \\
        \hat{a}_z
      \end{bmatrix}
      =
      {}^A\mathbf{C}^B
      \begin{bmatrix}
        \hat{b}_x \\
        \hat{b}_y \\
        \hat{b}_z
      \end{bmatrix}

:math:`\frac{{}^A\partial \bar{v}}{\partial q}`, ``v.diff(q, A)``
   Partial derivative of :math:`\bar{v}` with respect to :math:`q` when
   observed from :math:`A`.
:math:`\frac{{}^A d \bar{v}}{dt}`, ``v.dt(A)``
   Time derivative of :math:`\bar{v}` when observed from :math:`A`.
:math:`{}^A\bar{\omega}^B`, ``A_w_B``
   Angular velocity vector of reference frame or rigid body :math:`B` when
   observed from reference frame or rigid body :math:`A`.
:math:`{}^A\bar{\alpha}^B`, ``A_alp_B``
   Angular acceleration vector of reference frame or rigid body :math:`B` when
   observed from reference frame or rigid body :math:`A`.
:math:`\bar{r}^{P/O}`, ``r_O_P``
   Vector from point :math:`O` to point :math:`P`.
:math:`{}^A\bar{v}^P`, ``A_v_P``
   Translational velocity of point :math:`P` when observed from reference frame
   or rigid body :math:`A`.
:math:`{}^A\bar{a}^P`, ``A_a_P``
   Translational acceleration of point :math:`P` when observed from reference
   frame or rigid body :math:`A`.
:math:`\bar{f}_h(q_1, \ldots, q_n) = 0 \textrm{ where } \bar{f}_h \in \mathbb{R}^M`
   Vector function of :math:`M` holonomic constraint equations among the
   :math:`n` generalized coordinates.
