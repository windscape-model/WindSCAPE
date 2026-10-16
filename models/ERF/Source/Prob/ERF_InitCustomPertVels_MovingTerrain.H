    ParmParse pp_prob("prob");
    const MovingWave wave = MovingWave::from_inputs();
    const Real Ampl  = wave.a;
    const Real kp    = wave.k;
    const Real omega = wave.omega;

    // Optional start from a logarithmic wind u = (u*/kappa) ln((z + z0)/z0) over the surface, with a
    // small smooth 3-D seed that breaks the spanwise symmetry
    Real init_ustar    = zero;          pp_prob.query("init_ustar",    init_ustar);
    Real init_z0       = Real(2.0e-4);  pp_prob.query("init_z0",       init_z0);
    Real init_pert_amp = zero;          pp_prob.query("init_pert_amp", init_pert_amp);

    // 1: potential-flow field of the traveling wave (here and in the state initializer); 0: air at rest
    int  init_wave_field = 1;           pp_prob.query("init_wave_field", init_wave_field);
    const Real Ampl_field = (init_wave_field != 0) ? Ampl : zero;

    Real H           = geomdata.ProbHi()[2];
    Real Lx          = geomdata.ProbHi()[0] - geomdata.ProbLo()[0];
    Real Ly          = geomdata.ProbHi()[1] - geomdata.ProbLo()[1];
    Real l_kappa     = KAPPA;

    // x-velocity: potential flow of the wave, logarithmic wind and seed
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto prob_lo  = geomdata.ProbLo();
        const auto dx       = geomdata.CellSize();

        const Real x = prob_lo[0] + i * dx[0];
        const Real y = prob_lo[1] + (j + myhalf) * dx[1];
        Real z = fourth * (z_nd(i,j,k) + z_nd(i,j+1,k) + z_nd(i,j,k+1) + z_nd(i,j+1,k+1));

        Real z_base = Ampl * std::sin(kp * x);
        z -= z_base;

        Real fac     = std::cosh( kp * (z - H) ) / std::sinh(kp * H);
        Real u_wave  = -Ampl_field * omega * fac * std::sin(kp * x);

        Real u_mean  = (init_ustar > zero && z > zero) ?
                       (init_ustar / l_kappa) * std::log((z + init_z0) / init_z0) : zero;

        // Seed, zero at the surface and at the top
        Real fz   = (z > zero && z < H) ? std::sin(PI * z / H) : zero;
        Real u_pert = init_pert_amp * fz *
                      ( std::sin(two*PI*Real(3.0)*x/Lx) * std::cos(two*PI*Real(2.0)*y/Ly)
                      + std::sin(two*PI*Real(5.0)*x/Lx) * std::cos(two*PI*Real(4.0)*y/Ly) );

        x_vel_pert(i, j, k) = u_wave + u_mean + u_pert;
    });

    // y-velocity: seed only
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto prob_lo  = geomdata.ProbLo();
        const auto dx       = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + myhalf) * dx[0];
        const Real y = prob_lo[1] + j * dx[1];
        Real z = fourth * (z_nd(i,j,k) + z_nd(i+1,j,k) + z_nd(i,j,k+1) + z_nd(i+1,j,k+1));
        Real z_base = Ampl * std::sin(kp * x);
        z -= z_base;

        Real fz   = (z > zero && z < H) ? std::sin(PI * z / H) : zero;
        y_vel_pert(i, j, k) = init_pert_amp * fz *
                      ( std::cos(two*PI*Real(3.0)*x/Lx) * std::sin(two*PI*Real(2.0)*y/Ly)
                      + std::cos(two*PI*Real(5.0)*x/Lx) * std::sin(two*PI*Real(4.0)*y/Ly) );
    });

    // z-velocity: potential flow of the wave and seed
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto prob_lo  = geomdata.ProbLo();
        const auto dx         = geomdata.CellSize();

        Real x   = prob_lo[0] + (i + myhalf) * dx[0];
        Real y   = prob_lo[1] + (j + myhalf) * dx[1];
        Real z   = fourth * ( z_nd(i,j,k) + z_nd(i+1,j,k) + z_nd(i,j+1,k) + z_nd(i+1,j+1,k) );
        Real z_base = Ampl * std::sin(kp * x);

        z -= z_base;
        Real fac    = std::sinh( kp * (z - H) ) / std::sinh(kp * H);
        Real w_wave = Ampl_field * omega * fac * std::cos(kp * x);

        Real fz     = (z > zero && z < H) ? std::sin(PI * z / H) : zero;
        Real w_pert = init_pert_amp * Real(0.5) * fz *
                      std::sin(two*PI*Real(4.0)*x/Lx) * std::sin(two*PI*Real(3.0)*y/Ly);

        z_vel_pert(i, j, k) = w_wave + w_pert;
    });
