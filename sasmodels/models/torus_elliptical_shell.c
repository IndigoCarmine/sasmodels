static double form_volume(double radius, double core_radius, double thickness,
                          double nu_core, double nu_shell) {
  double ao = core_radius + thickness;
  double bo = nu_core * core_radius + nu_shell * thickness;
  double area = M_PI * ao * bo;
  return 2.0 * M_PI * radius * area;
}

static double radius_from_volume(double core_radius, double thickness,
                                 double nu_core, double nu_shell) {
  return cbrt(form_volume(1.0, core_radius, thickness, nu_core, nu_shell) /
              M_4PI_3);
}

static double radius_from_diagonal(double radius, double core_radius,
                                   double thickness, double nu_core,
                                   double nu_shell) {
  return radius + core_radius + thickness;
}

static double max_radius(double core_radius, double thickness, double nu_core,
                         double nu_shell, double torus_radius) {
  // equatorial radius of the outer ellipse
  double r_e = core_radius + thickness;

  // polar radius of the outer ellipse
  double r_p = nu_core * core_radius + nu_shell * thickness;

  double denom = square(r_p) - square(r_e);

  if (denom > 0.0) {
    double cos_theta = (torus_radius * r_e) / denom;

    if (fabs(cos_theta) <= 1.0) {
      double nu = r_p / r_e;

      double r_max_sq = square(torus_radius) + square(r_e) +
                        square(r_p) * ((square(nu) - 1.0) / square(nu));

      return sqrt(r_max_sq);
    }
  }

  return torus_radius + r_e;
}

static double radius_effective(int mode, double radius, double core_radius,
                               double thickness, double nu_core,
                               double nu_shell) {
  switch (mode) {
    case 1:
      return radius_from_diagonal(radius, core_radius, thickness, nu_core,
                                  nu_shell);
    case 2:
      return radius_from_volume(core_radius, thickness, nu_core, nu_shell);
    case 3:
      return max_radius(core_radius, thickness, nu_core, nu_shell, radius);
    case 4:
    default:
      return radius;
  }
}

static double F_torus(double Q, double theta, double R, double x, double nu,
                      double delta_eta) {
  // integral_{R - x}^{R + x}
  // 4 * pi * r * J0(Q * r * sin(theta)) * ganmma
  // sinx_x(Q * ganmma * cos(theta)) dr

  // Set lower integration bound to R-x (not min(R-x,0)), since a torus has a
  // hole and R-x > 0 is always valid.

  double gamma = 0, int_r_delta = 0, r = 0, f_total = 0;
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

static void Fq(double q, double* F1, double* F2, double radius,
               double core_radius, double thickness, double nu_core,
               double nu_shell, double sld_core, double sld_shell,
               double sld_solvent) {
  // F2 = integral_{0}^{pi / 2}
  //  | F_torus(Q, theta, R, core_radius + thickness, nu, sld_shell -
  //  sld_solvent)   // shell
  //    - F_torus(Q, theta, R, core_radius, nu, sld_core - sld_solvent) //
  //    core
  //  |^2 dtheta
  double F_diff = 0.0, theta = 0.0, F1_total = 0.0, F2_total = 0.0;

  double nu_outer =
      (nu_shell * thickness + nu_core * core_radius) /
      (thickness + core_radius);  // aspect ratio of the outer ellipse

  for (int i = 0; i < GAUSS_N; i++) {
    // translate a point in[-1, 1] to a point in[0, pi / 2]
    theta = GAUSS_Z[i] * M_PI_4 + M_PI_4;

    F_diff =
        F_torus(q, theta, radius, core_radius + thickness, nu_outer,
                sld_shell - sld_solvent) -
        F_torus(q, theta, radius, core_radius, nu_core, sld_shell - sld_core);

    F1_total += GAUSS_W[i] * F_diff * sin(theta);
    F2_total += GAUSS_W[i] * square(F_diff) * sin(theta);
  }

  // multiply by pi/4 to get the integral over [0, pi/2]

  *F1 = 1e-2 * F1_total * M_PI_4;
  *F2 = 1e-4 * F2_total * M_PI_4;
}

static double Iq(double q, double radius, double core_radius, double thickness,
                 double nu_core, double nu_shell, double sld_core,
                 double sld_shell, double sld_solvent) {
  double F1 = 0.0, F2 = 0.0;
  Fq(q, &F1, &F2, radius, core_radius, thickness, nu_core, nu_shell, sld_core,
     sld_shell, sld_solvent);
  return F2;
}