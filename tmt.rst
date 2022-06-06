=====================================================
Unconstrained Equatiosn of Motion With the TMT Method
=====================================================

There are many mathematical methods available to formulate the equations of
motion given a description of a multibody system model. The different methods
offer different advantages and disadvantages over the first methods proposed by
Netwon and Euler. For example, Joseph-Louis Lagrange developed a way to arrive
at the equations of motion from the descriptions of kinetic and potential
energy of the system. Later Sir William Hamilton reformulated Lagrange's
approach in terms of generalized momenta. In this chapter, we present an
alternative method presented in [Vallery2020]_ call the "TMT Method".

Vallery and Schwab present a matrix that can be used to directly transform the
mass and inertia of a system of rigid bodies expressed in a Cartesian
coordinate system into the mass matrix for a specific set of generalized
coordinates. This matrix is populated by the measure numbers of the partial
velocities expressed in the inertial reference frame. The method has some
advantages:

- Only the position and orientation of the bodies need be determined.
- Contributing forces can be formaluated as generalized forces, or not, or
  both.



Given :math:`\nu` rigid bodies in a multibody system, the velocities of each
mass center and the angular velocities of each body in an inertial reference
frame :math:`N` can be written in column vector :math:`\bar{v}` form by
extracting the measure numbers of each velocity term.

.. math::

   \bar{v}(\dot{\bar{q}}, \bar{q}, t) =
   \begin{bmatrix}
   {}^N\bar{v}^{B_{1o}} \cdot \hat{n}_x \\
   {}^N\bar{v}^{B_{1o}} \cdot \hat{n}_y \\
   {}^N\bar{v}^{B_{1o}} \cdot \hat{n}_z \\
   {}^N\bar{\omega}^{B_1} \cdot \hat{n}_x \\
   {}^N\bar{\omega}^{B_1} \cdot \hat{n}_y \\
   {}^N\bar{\omega}^{B_1} \cdot \hat{n}_z \\
   \vdots \\
   {}^N\bar{v}^{B_{\nu o}} \cdot \hat{n}_x \\
   {}^N\bar{v}^{B_{\nu o}} \cdot \hat{n}_y \\
   {}^N\bar{v}^{B_{\nu o}} \cdot \hat{n}_z \\
   {}^N\bar{\omega}^{B_\nu} \cdot \hat{n}_x \\
   {}^N\bar{\omega}^{B_\nu} \cdot \hat{n}_y \\
   {}^N\bar{\omega}^{B_\nu} \cdot \hat{n}_z \\
   \end{bmatrix}
   \in
   \mathbb{R}^{6\nu}

The measure numbers of the partial velocities with respect to each speed in
:math:`\dot{\bar{q}}` can be efficiently found with the Jacobian and stored in
a matrix :math:`\mathbf{T}`.

.. math::

   \mathbf{T} = J_\bar{v},\dot{\bar{q}}

The mass an inertia of each body can be expressed in a 

.. math::

   \mathbf{M}_{B_1} =
   \begin{bmatrix}
   m_{B_1} & 0 & 0 & 0 & 0 & 0 \\
   0 & m_{B_1} & 0 & 0 & 0 & 0 \\
   0 & 0 & m_{B_1} & 0 & 0 & 0 \\
   0 & 0 & 0 &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_x\hat{n}_x &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_x\hat{n}_y &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_x\hat{n}_z \\
   0 & 0 & 0 &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_y\hat{n}_x &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_y\hat{n}_y &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_y\hat{n}_z \\
   0 & 0 & 0 &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_z\hat{n}_x &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_z\hat{n}_y &
   \breve{I}^{B_1/B_{1o}} \cdot \hat{n}_z\hat{n}_z \\
   \end{bmatrix}

.. math::

   \mathbf{M} =
   \begin{bmatrix}
   \mathbf{M}_{B_1} & \mathbf{0}       & \ldots     & \mathbf{0} \\
   \mathbf{0}       & \mathbf{M}_{B_2} & \ldots     & \vdots \\
   \vdots           & \vdots           & \ddots     & \vdots \\
   \mathbf{0}       & \mathbf{0}       & \ldots     & \mathbf{M}_{B_\nu}
   \end{bmatrix}

Vallery and Schwab show that the the mass matrix of the system is:

.. math::

   \mathbf{M}_d = \mathbf{T}^T \mathbf{M} \mathbf{T}

Now given a vector :math:`\bar{F}` of resultant forces and torques of couples acting on each rigid body:

.. math::

   \bar{F} =
   \begin{bmatrix}
   \bar{R}^{B_{1o}} \cdot \hat{n}_x \\
   \bar{R}^{B_{1o}} \cdot \hat{n}_y \\
   \bar{R}^{B_{1o}} \cdot \hat{n}_z \\
   \bar{T}^{B_1} \cdot \hat{n}_x \\
   \bar{T}^{B_1} \cdot \hat{n}_y \\
   \bar{T}^{B_1} \cdot \hat{n}_z \\
   \vdots \\
   \bar{R}^{B_{2o}} \cdot \hat{n}_x \\
   \bar{R}^{B_{2o}} \cdot \hat{n}_y \\
   \bar{R}^{B_{2o}} \cdot \hat{n}_z \\
   \bar{T}^{B_2} \cdot \hat{n}_x \\
   \bar{T}^{B_2} \cdot \hat{n}_y \\
   \bar{T}^{B_2} \cdot \hat{n}_z \\
   \end{bmatrix}

.. math::

   \bar{g} = \frac{d\bar{v}}{dt}\bigg\rvert_{\ddot{\bar{q}}=\bar{0}}

.. math::

   \mathbf{T}^T \mathbf{M} \mathbf{T} \ddot{\bar{q}} =
   \mathbf{T}^T\left(\bar{F} - \mathbf{M}\bar{g}\right)


