#ifdef USE_OPENCL
__constant double A =
    109.5 / 2.0 * M_PI / 180.0;  // half of the angle between arms in radians
#else
const double A =
    109.5 / 2.0 * M_PI / 180;  // half of the angle between arms in radians
#endif

static double u_n(int n, double theta, double alpha) {
  const double phi[4] = {0.0, M_PI_2, M_PI, 3.0 * M_PI_2};
  const double sign[4] = {1.0, -1.0, 1.0, -1.0};
  return sign[n] * cos(A) * cos(theta) +
         sin(A) * sin(theta) * cos(alpha - phi[n]);
}

static double Fq_n(double q, double u, double L, double R) {
  double quL2 = q * u * L * 0.5;

  double mu = sqrt(fmax(0.0, 1.0 - u * u));
  double qmuR = q * mu * R;

  return sas_sinx_x(quL2) * sas_2J1x_x(qmuR);
}

static double form_volume(double L, double R) {
  // V = 4 * \pi * R^2 * L
  return 4.0 * M_PI * R * R * L;
}
static double Iq_same(double q, double L, double R, double delta_sld) {
  double total = 0.0;

  for (int dtheta = 0; dtheta < GAUSS_N; dtheta++) {
    double theta =
        M_PI_2 * (GAUSS_Z[dtheta] + 1.0);  // map from [-1, 1] to [0, pi]
    double w_theta =
        GAUSS_W[dtheta] * M_PI_2;  // adjust weight for the new range

    double integral_alpha = 0.0;
    for (int dalpha = 0; dalpha < GAUSS_N; dalpha++) {
      double alpha =
          M_PI * (GAUSS_Z[dalpha] + 1.0);  // map from [-1, 1] to [0, 2*pi]
      double w_alpha =
          GAUSS_W[dalpha] * M_PI;  // adjust weight for the new range

      double sum_arms = 0.0;
      for (int n = 0; n < 4; n++) {
        double u = u_n(n, theta, alpha);

        double fn = Fq_n(q, u, L, R);
        sum_arms += fn * fn;
      }
      integral_alpha += sum_arms * w_alpha;
    }
    total += integral_alpha * w_theta * sin(theta);
  }
  return delta_sld * delta_sld * total / M_PI / M_PI;
}

static double Iq_different(double q, double L, double R, double delta_sld) {
  double total = 0.0;
  const int nm_pairs[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

  for (int dtheta = 0; dtheta < GAUSS_N; dtheta++) {
    double theta =
        M_PI_2 * (GAUSS_Z[dtheta] + 1.0);  // map from [-1, 1] to [0, pi]
    double w_theta =
        GAUSS_W[dtheta] * M_PI_2;  // adjust weight for the new range

    double integral_alpha = 0.0;
    for (int dalpha = 0; dalpha < GAUSS_N; dalpha++) {
      double alpha =
          M_PI * (GAUSS_Z[dalpha] + 1.0);  // map from [-1, 1] to [0, 2*pi]
      double w_alpha =
          GAUSS_W[dalpha] * M_PI;  // adjust weight for the new range

      double sum_arms = 0.0;
      for (int i = 0; i < 6; i++) {
        int n = nm_pairs[i][0];
        int m = nm_pairs[i][1];
        double u = u_n(n, theta, alpha);
        double v = u_n(m, theta, alpha);
        sum_arms +=
            Fq_n(q, u, L, R) * Fq_n(q, v, L, R) * cos(q * (u - v) * L / 2.0);
      }
      integral_alpha += sum_arms * w_alpha;
    }
    total += integral_alpha * w_theta * sin(theta);
  }
  return delta_sld * delta_sld * total * 3 / M_PI / M_PI;
}

static double Iq(double q, double L, double R, double sld_particle,
                 double sld_solvent) {
  double contrast = sld_particle - sld_solvent;
  return Iq_same(q, L, R, contrast) + Iq_different(q, L, R, contrast);
}

// static double Iq_old(double q, double L, double R, double sld_particle,
//                      double sld_solvent) {
//   double contrast = sld_particle - sld_solvent;
//   double total = 0.0;

//   for (int dtheta = 0; dtheta < GAUSS_N; dtheta++) {
//     double theta =
//         M_PI_2 * (GAUSS_Z[dtheta] + 1.0);  // map from [-1, 1] to [0, pi]
//     double w_theta =
//         GAUSS_W[dtheta] * M_PI_2;  // adjust weight for the new range

//     double integral_alpha = 0.0;
//     for (int dalpha = 0; dalpha < GAUSS_N; dalpha++) {
//       double alpha =
//           M_PI * (GAUSS_Z[dalpha] + 1.0);  // map from [-1, 1] to [0, 2*pi]
//       double w_alpha =
//           GAUSS_W[dalpha] * M_PI;  // adjust weight for the new range

//       double sum_arms = 0.0;
//       for (int n = 0; n < 4; n++) {
//         for (int m = 0; m < 4; m++) {
//           double u = u_n(n, theta, alpha);
//           double v = u_n(m, theta, alpha);
//           sum_arms +=
//               Fq_n(q, u, L, R) * Fq_n(q, v, L, R) * cos(q * (u - v) * L
//               / 2.0);
//         }
//       }
//       integral_alpha += sum_arms * w_alpha * sin(theta);
//     }
//     total += integral_alpha * w_theta;
//   }
//   return contrast * contrast * total * 2 / M_PI;
// }