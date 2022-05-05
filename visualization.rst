===============================
Three Dimensional Visualization
===============================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`visualization` or Jupyter Notebook:
   :jupyter-download:notebook:`visualization`.

.. jupyter-execute::

   from scipy.integrate import solve_ivp
   import numpy as np
   import sympy as sm
   import sympy.physics.mechanics as me

.. admonition:: Modeling and Simulation Code
   :class: dropdown

   .. jupyter-execute::

      m, g, kt, kl, l = sm.symbols('m, g, k_t, k_l, l')
      q1, q2, q3 = me.dynamicsymbols('q1, q2, q3')
      u1, u2, u3 = me.dynamicsymbols('u1, u2, u3')

      N = me.ReferenceFrame('N')
      A = me.ReferenceFrame('A')
      B = me.ReferenceFrame('B')

      A.orient_axis(N, q1, N.z)
      B.orient_axis(A, q2, A.x)

      A.set_ang_vel(N, u1*N.z)
      B.set_ang_vel(A, u2*A.x)

      O = me.Point('O')
      Ao = me.Point('A_O')
      Bo = me.Point('B_O')
      Q = me.Point('Q')

      Ao.set_pos(O, l/2*A.x)
      Bo.set_pos(O, l*A.x)
      Q.set_pos(Bo, q3*B.y)

      O.set_vel(N, 0)
      Ao.v2pt_theory(O, N, A)
      Bo.v2pt_theory(O, N, A)
      Q.set_vel(B, u3*B.y)
      Q.v1pt_theory(Bo, N, B)

      t = me.dynamicsymbols._t

      qdot_repl = {q1.diff(t): u1,
                   q2.diff(t): u2,
                   q3.diff(t): u3}

      Q.set_acc(N, Q.acc(N).xreplace(qdot_repl))

      R_Ao = m*g*N.x
      R_Bo = m*g*N.x + kl*q3*B.y
      R_Q = m/4*g*N.x - kl*q3*B.y
      T_A = -kt*q1*N.z + kt*q2*A.x
      T_B = -kt*q2*A.x

      I = m*l**2/12
      I_A_Ao = I*me.outer(A.y, A.y) + I*me.outer(A.z, A.z)
      I_B_Bo = I*me.outer(B.x, B.x) + I*me.outer(B.z, B.z)

      points = [Ao, Bo, Q]
      forces = [R_Ao, R_Bo, R_Q]
      masses = [m, m, m/4]

      frames = [A, B]
      torques = [T_A, T_B]
      inertias = [I_A_Ao, I_B_Bo]

      Fr_bar = []
      Frs_bar = []

      for ur in [u1, u2, u3]:

         Fr = 0
         Frs = 0

         for Pi, Ri, mi in zip(points, forces, masses):
            vr = Pi.vel(N).diff(ur, N)
            Fr += vr.dot(Ri)
            Rs = -mi*Pi.acc(N)
            Frs += vr.dot(Rs)

         for Bi, Ti, Ii in zip(frames, torques, inertias):
            wr = Bi.ang_vel_in(N).diff(ur, N)
            Fr += wr.dot(Ti)
            Ts = -(Bi.ang_acc_in(N).dot(Ii) +
                   me.cross(Bi.ang_vel_in(N), Ii).dot(Bi.ang_vel_in(N)))
            Frs += wr.dot(Ts)

         Fr_bar.append(Fr)
         Frs_bar.append(Frs)

      Fr = sm.Matrix(Fr_bar)
      Frs = sm.Matrix(Frs_bar)

      q = sm.Matrix([q1, q2, q3])
      u = sm.Matrix([u1, u2, u3])
      p = sm.Matrix([g, kl, kt, l, m])

      qd = q.diff(t)
      ud = u.diff(t)

      ud_zerod = {udr: 0 for udr in ud}

      Mk = -sm.eye(3)
      gk = u

      Md = Frs.jacobian(ud)
      gd = Frs.xreplace(ud_zerod) + Fr

      eval_eom = sm.lambdify((q, u, p), [Mk, gk, Md, gd])

      def eval_rhs(t, x, p):
          """Return the right hand side of the explicit ordinary differential
          equations which evaluates the time derivative of the state ``x`` at time
          ``t``.

          Parameters
          ==========
          t : float
             Time in seconds.
          x : array_like, shape(6,)
             State at time t: [q1, q2, q3, u1, u2, u3]
          p : array_like, shape(5,)
             Constant parameters: [g, kl, kt, l, m]

          Returns
          =======
          xd : ndarray, shape(6,)
              Derivative of the state with respect to time at time ``t``.

          """

          # unpack the q and u vectors from x
          q = x[:3]
          u = x[3:]

          # evaluate the equations of motion matrices with the values of q, u, p
          Mk, gk, Md, gd = eval_eom(q, u, p)

          # solve for q' and u'
          qd = np.linalg.solve(-Mk, np.squeeze(gk))
          ud = np.linalg.solve(-Md, np.squeeze(gd))

          # pack dq/dt and du/dt into a new state time derivative vector dx/dt
          xd = np.empty_like(x)
          xd[:3] = qd
          xd[3:] = ud

          return xd

      q_vals = np.array([
          np.deg2rad(25.0),  # q1, rad
          np.deg2rad(5.0),  # q2, rad
          0.1,  # q3, m
      ])

      u_vals = np.array([
          0.1,  # u1, rad/s
          2.2,  # u2, rad/s
          0.3,  # u3, m/s
      ])

      p_vals = np.array([
          9.81,  # g, m/s**2
          2.0,  # kl, N/m
          0.01,  # kt, Nm/rad
          0.6,  # l, m
          1.0,  # m, kg
      ])

      x0 = np.empty(6)
      x0[:3] = q_vals
      x0[3:] = u_vals

      fps = 50
      t0, tf = 0.0, 10.0
      ts = np.linspace(t0, tf, num=int(fps*(tf - t0)))
      result = solve_ivp(eval_rhs, (t0, tf), x0, args=(p_vals,), t_eval=ts)
      xs = result.y.T

.. jupyter-execute::

   ts

.. jupyter-execute::

   xs

.. jupyter-execute::

   import pythreejs as p3js

   p3js.CylinderBufferGeometry(
       radiusTop=5,
       radiusBottom=10,
       height=50,
       radialSegments=20,
       heightSegments=10,
       openEnded=False,
       thetaStart=0,
       thetaLength=2.0*np.pi)

Lot's of geometry types https://pythreejs.readthedocs.io/en/stable/examples/Geometries.html

https://en.wikipedia.org/wiki/Transformation_matrix

.. math::

   \mathbf{T} = \begin{bmatrix}
   {}^N\mathbf{C}^B & \bar{0} \\
   \bar{r}^{B_o/O} & 0
   \end{bmatrix}

.. jupyter-execute::

   TA = sm.eye(4)
   A_ = me.ReferenceFrame('A_')
   A_.orient_axis(A, sm.pi/2, A.z)
   TA[:3, :3] = A_.dcm(N)
   TA[3, :3] = sm.transpose(Ao.pos_from(O).to_matrix(N))

   TA

.. jupyter-execute::

   TB = sm.eye(4)

   TB[:3, :3] = B.dcm(N)
   TB[3, :3] = sm.transpose(Bo.pos_from(O).to_matrix(N))

   TB

.. jupyter-execute::

   TQ = sm.eye(4)

   TQ[:3, :3] = B.dcm(N)
   TQ[3, :3] = sm.transpose(Q.pos_from(O).to_matrix(N))

   TQ

.. jupyter-execute::

   eval_transform = sm.lambdify((q, p), (TA, TB, TQ))
   eval_transform(q_vals, p_vals)

.. jupyter-execute::

   transforms_A = []
   transforms_B = []
   transforms_Q = []
   for qi in result.y.T:
      one, two, three = eval_transform(qi[:3], p_vals)
      transforms_A.append(one.flatten())
      transforms_B.append(two.flatten())
      transforms_Q.append(three.flatten())

   TA_vals = np.array(transforms_A).tolist()
   TB_vals = np.array(transforms_B).tolist()
   TQ_vals = np.array(transforms_Q).tolist()

.. jupyter-execute::

   cylA_geom = p3js.CylinderBufferGeometry(
       radiusTop=p_vals[3]/20,
       radiusBottom=p_vals[3]/20,
       height=p_vals[3],  # l
       #radialSegments=20,
       )
   cylA_geom

.. jupyter-execute::

   cylB_geom = p3js.CylinderBufferGeometry(
       radiusTop=p_vals[3]/20,
       radiusBottom=p_vals[3]/20,
       height=p_vals[3],  # l
       #radialSegments=20,
       )
   cylB_geom

.. jupyter-execute::

   sphQ_geom = p3js.SphereBufferGeometry(
        radius=p_vals[3]/16)
   sphQ_geom

.. jupyter-execute::

   arrow_length = 0.3

   rodA = p3js.Mesh(
       geometry=cylA_geom,
       material=p3js.MeshStandardMaterial(color='red'),
       name='rodA',
   )
   rodA.matrixAutoUpdate = False
   rodA.add(p3js.AxesHelper(arrow_length))
   rodA.matrix = TA_vals[0]

   rodB = p3js.Mesh(
       geometry=cylB_geom,
       material=p3js.MeshStandardMaterial(color='blue'),
       name='rodB',
   )
   rodB.matrixAutoUpdate = False
   rodB.add(p3js.AxesHelper(arrow_length))
   rodB.matrix = TB_vals[0]

   pointQ = p3js.Mesh(
       geometry=sphQ_geom,
       material=p3js.MeshStandardMaterial(color='green'),
       name='pointQ',
   )
   pointQ.matrixAutoUpdate = False
   pointQ.add(p3js.AxesHelper(arrow_length))
   pointQ.matrix = TQ_vals[0]

The X axis is red. The Y axis is green. The Z axis is blue.

.. jupyter-execute::

   view_width = 600
   view_height = 400

   camera = p3js.PerspectiveCamera(position=[1, 0.6, 1],
                                  aspect=view_width/view_height)
   camera.up = (-1, 0, 0)
   key_light = p3js.DirectionalLight(position=[0, 10, 10])

   transform_track_rodA = p3js.VectorKeyframeTrack(
       name="scene/rodA.matrix",
       times=result.t,
       values=TA_vals
   )

   transform_track_rodB = p3js.VectorKeyframeTrack(
       name="scene/rodB.matrix",
       times=result.t,
       values=TB_vals
   )

   transform_track_pointQ = p3js.VectorKeyframeTrack(
       name="scene/pointQ.matrix",
       times=result.t,
       values=TQ_vals
   )

   ambient_light = p3js.AmbientLight()
   camera_clip = p3js.AnimationClip(tracks=[transform_track_rodB,
   transform_track_rodA, transform_track_pointQ], duration=result.t[-1])

   axes = p3js.AxesHelper()
   scene = p3js.Scene(children=[rodA, rodB, pointQ, axes, camera, key_light, ambient_light])
   controller = p3js.OrbitControls(controlling=camera)
   renderer = p3js.Renderer(camera=camera, scene=scene, controls=[controller],
                       width=view_width, height=view_height)

   camera_action = p3js.AnimationAction(p3js.AnimationMixer(scene), camera_clip, scene)

https://pythreejs.readthedocs.io/en/stable/examples/Animation.html

.. jupyter-execute::

   renderer

.. jupyter-execute::

   camera_action

