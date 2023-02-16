========
Notation
========

This explains the notation in the book. The mathematical symbol is shown and
then an example of a variable name that we use in the code.

.. todo:: Add the notation that is in the hand drawn figures.

:math:`x`, ``x``, |notation-scalar|
   Scalars are normal mathematical font.
:math:`\bar{v}`, ``v``
   Vectors are indicated with a bar.
:math:`\hat{u}`, ``uhat = u.normalize()``, |notation-unit-vec|
   Unit vectors are indicated with a hat.
:math:`A`, ``A``, |notation-ref-frame|
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
:math:`{}^A\mathbf{C}^B`, ``A_C_B``, |notation-dcm|
   Direction cosine matrix relating reference frames (or rigid bodies)
   :math:`B` and :math:`A` where this relation between the right handed
   mutually perpendicular unit vectors fixed in the two reference frames follow
   this relationship:

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

:math:`\frac{{}^A\partial \bar{v}}{\partial q}`, ``dvdqA = v.diff(q, A)``
   Partial derivative of :math:`\bar{v}` with respect to :math:`q` when
   observed from :math:`A`.
:math:`\frac{{}^A d \bar{v}}{dt}`, ``dvdtA = v.dt(A)``
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
:math:`\bar{f}_h(q_1, \ldots, q_N, t) = 0 \textrm{ where } \bar{f}_h \in \mathbb{R}^M`, ``fh``
   Vector function of :math:`M` holonomic constraint equations among the
   :math:`N` coordinates.
:math:`\bar{f}_n(u_1, \ldots, u_n, q_1, \ldots, q_n, t) = 0 \textrm{ where } \bar{f}_n \in \mathbb{R}^m`, ``fn``
   Vector function of :math:`m` nonholonomic constraint equations among the
   :math:`n` generalized speeds and generalized coordinates.
:math:`\bar{I}^{B/O}_a`, ``I_B_O_a``
   Inertia vector of rigid body :math:`B` with respect to point :math:`O` about
   the unit vector :math:`\hat{n}_a`.
:math:`\breve{Q}`, ``Q``
   Dyadics are indicated with a breve accent.
:math:`\breve{I}^{B/O}`, ``I_B_O``
   Inertia dyadic of body :math:`B` or set of particles :math:`B` with respect
   to point :math:`O`.
:math:`\breve{I}^{B/B_o}`, ``I_B_Bo``
   Central inertia dyadic of body :math:`B` or set of particles :math:`B` with respect
   to mass center :math:`B_o`.
:math:`{}^A \bar{H}^{B/O}`, ``A_H_B_O``
   Angular momentum of rigid body :math:`B` with respect to point :math:`O` in
   reference frame :math:`A`.
:math:`\bar{R}^{S}`, ``R_S``
   Resultant of the vector set :math:`S`.
:math:`\bar{R}^{S/Q}`, ``R_S_Q``
   Resultant of the vector set :math:`S` bound to a line of action through
   point :math:`Q`.
:math:`\bar{M}^{S/P}`, ``M_S_P``
   Moment of the resultant of the vector set :math:`S` about point :math:`P`.
:math:`\bar{T}^{B}`, ``T_B``
   Torque of couple acting on reference frame or body :math:`B`.
:math:`{}^A\bar{v}_r^P`, ``v_P_r``
   r\ :sup:`th` holonomic partial velocity of point :math:`P` in reference
   frame :math:`A` associated with the generalized speed :math:`u_r`.
:math:`{}^A\bar{\omega}_r^B`, ``w_B_r``
   r\ :sup:`th` holonomic partial angular velocity of reference frame :math:`B`
   in reference frame :math:`A` associated with the generalized speed
   :math:`u_r`.
:math:`{}^A\tilde{v}_r^P`, ``v_P_r``
   r\ :sup:`th` nonholonomic partial velocity of point :math:`P` in reference
   frame :math:`A` associated with the generalized speed :math:`u_r`.
:math:`{}^A\tilde{\omega}_r^B`, ``w_B_r``
   r\ :sup:`th` nonholonomic partial angular velocity of reference frame
   :math:`B` in reference frame :math:`A` associated with the generalized speed
   :math:`u_r`.
:math:`F_r`, ``F1``
   r\ :sup:`th` holonomic generalized active force associated with the
   generalized speed :math:`u_r`.
:math:`\tilde{F}_r`, ``F1``
   r\ :sup:`th` nonholonomic generalized active force associated with the
   generalized speed :math:`u_r`.
:math:`\bar{F}_r`, ``Fr``
   Column vector of all generalized active forces (holonomic or nonholonomic).
:math:`F^*_r`, ``F1s``
   r\ :sup:`th` holonomic generalized inertia force associated with the
   generalized speed :math:`u_r`.
:math:`\tilde{F}^*_r`, ``F1s``
   r\ :sup:`th` nonholonomic generalized inertia force associated with the
   generalized speed :math:`u_r`.
:math:`\bar{F}^*_r`, ``Frs``
   Column vector of all generalized active forces (holonomic or nonholonomic).
:math:`\mathbf{J}_{\bar{v},\bar{u}}`
   The Jacobian of the vector function :math:`\bar{v}` with respect to the
   entries in vector :math:`\bar{u}` where the :math:`(i,j)` entries of the
   Jacobian are :math:`\mathbf{J}_{ij} = \frac{\partial v_i}{\partial u_j}`.

.. |notation-scalar| image:: figures/notation-scalar.svg
   :height: 10px

.. |notation-unit-vec| image:: figures/notation-unit-vec.svg
   :height: 15px

.. |notation-ref-frame| image:: figures/notation-ref-frame.svg
   :height: 20px

.. |notation-dcm| image:: figures/notation-dcm.svg
   :height: 20px
