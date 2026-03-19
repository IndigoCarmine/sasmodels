static double form_volume(double radius, double core_radius, double thickness,
                          double nu) {
  double a = core_radius;
  double b = nu * a;
  double ao = a + thickness;
  double bo = nu * ao;
  double area = M_PI * ao * bo - M_PI * a * b;
  return 2.0 * M_PI * radius * area;
}

static double F_torus(double Q, double theta, double R, double x, double nu,
                      double delta_eta) {
  // integral_{R - x}^{R + x}
  // 4 * pi * r * J0(Q * r * sin(theta)) * ganmma
  // sinx_x(Q * ganmma * cos(theta)) dr

  // Set lower integration bound to R-x (not min(R-x,0)), since a torus has a
  // hole and R-x > 0 is always valid.

  double gamma, int_r_delta, r, f_total, gamma_arg;
  const double square_x = square(x);
  const double Q_sin_theta = Q * sin(theta);
  const double Q_cos_theta = Q * cos(theta);

  f_total = 0.0;

  for (int i = 0; i < GAUSS_N; i++) {
    // translate a point in[-1, 1] to a point in[R - x, R + x]
    r = GAUSS_Z[i] * x + R;

    gamma = nu * sqrt(square_x - square(r - R));

    int_r_delta =
        r * sas_J0(r * Q_sin_theta) * sas_sinx_x(Q_cos_theta * gamma) * gamma;
    f_total += GAUSS_W[i] * int_r_delta;
  }

  return 4.0 * M_PI * f_total * x *
         delta_eta;  // multiply by x to get the integral over [R - x, R + x]
}

static double Iq(double q, double radius, double core_radius, double thickness,
                 double nu, double sld_core, double sld_shell,
                 double sld_solvent) {
  // integral_{0}^{pi / 2}
  //  | F_torus(Q, theta, R, core_radius + thickness, nu, sld_shell -
  //  sld_solvent)   // shell
  //    - F_torus(Q, theta, R, core_radius, nu, sld_core - sld_solvent) //
  //    core
  //  |^2 dtheta
  double F_diff, theta;
  double I_total = 0.0;

  for (int i = 0; i < GAUSS_N; i++) {
    // translate a point in[-1, 1] to a point in[0, pi / 2]
    theta = GAUSS_Z[i] * M_PI_4 + M_PI_4;

    F_diff = F_torus(q, theta, radius, core_radius + thickness, nu,
                     sld_shell - sld_solvent) -
             F_torus(q, theta, radius, core_radius, nu, sld_shell - sld_core);

    I_total += GAUSS_W[i] * square(F_diff) * sin(theta);
  }

  return I_total *
         M_PI_4;  // multiply by pi/4 to get the integral over [0, pi/2]
}