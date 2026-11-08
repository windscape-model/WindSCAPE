
    // Parse params
    ParmParse pp("prob");
    const MovingWave wave = MovingWave::from_inputs();
    const Real Ampl  = wave.a;
    const Real kp    = wave.k;
    const Real omega = wave.omega;

    // 1: pressure of the potential flow of the wave; 0: air at rest
    int  init_wave_field = 1; pp.query("init_wave_field", init_wave_field);
    const Real Ampl_field = (init_wave_field != 0) ? Ampl : zero;

    Real H           = geomdata.ProbHi()[2];

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto prob_lo  = geomdata.ProbLo();
        const auto dx       = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + myhalf) * dx[0];

        Real z = z_cc(i,j,k);
        Real z_base = Ampl * std::sin(kp * x);
        z -= z_base;

        Real fac     = std::cosh( kp * (z - H) ) / std::sinh(kp * H);
        Real p_prime = -(Ampl_field * omega * omega / kp) * fac * std::sin(kp * x) * r_hse(i,j,k);
        Real p_total = p_hse(i,j,k) + p_prime;

        // Define (rho theta) given pprime
        state_pert(i, j, k, RhoTheta_comp) = getRhoThetagivenP(p_total) - getRhoThetagivenP(p_hse(i,j,k));

        // Set scalar = 0 everywhere
        state_pert(i, j, k, RhoScalar_comp) = state_pert(i,j,k,Rho_comp);
    });
