void productionRate(
    double* wdot,
    const double* sc,
    const double T,
    const double Te,
    double EN,
    double* enerExch)
{

    double tc[5] = {
        log(T), T, T * T, T * T * T, T * T * T * T}; // temperature cache
    const double invT = 1.0 / tc[1];
    const double logT = log(T / 300.0);

    // reference concentration: P_atm / (RT) in inverse mol/m^3
    const double refC = 101325 / 8.31446 * invT;
    const double refCinv = 1 / refC;

    for (int i = 0; i < 44; ++i)
    {
        wdot[i] = 0.0;
    }

    // compute the mixture concentration
    double mixture = 0.0;
    for (int i = 0; i < 44; ++i)
    {
        mixture += sc[i];
    }

    // compute the Gibbs free energy
    double g_RT[44];
    gibbs(g_RT, T);

    // Precalculating values for electron energy exchange evaluation
    double ne = sc[E_ID] * 6.02214085774e23;
    double Ue = 1.5 * Te * ne * 1.380649e-23;

    // NOTE: units of JANEV fits for electron impact rxns are cm3/s and must be
    // convered to m3/mol-s Precalculating values
    double Janev_sum;
    double invTe = (Te == 0) ? 1.0 : 1.0 / Te;
    double TeeV = Te / 11595.0;
    double logTe =
    log(TeeV); // Fits are performed assuming Te is eV rather than K
    double invTeeV = (Te == 0) ? 1.0 : 1.0 / (TeeV);
    double Te_pow[] = {
        1.0,
        logTe,
        amrex::Math::powi<2>(logTe),
        amrex::Math::powi<3>(logTe),
        amrex::Math::powi<4>(logTe),
        amrex::Math::powi<5>(logTe),
        amrex::Math::powi<6>(logTe),
        amrex::Math::powi<7>(logTe),
        amrex::Math::powi<8>(logTe)};
    double invTe_pow[] = {
        invTeeV, amrex::Math::powi<2>(invTeeV), amrex::Math::powi<3>(invTeeV),
        amrex::Math::powi<4>(invTeeV)};

    tc[0] = log(T / 300.0); // temperature factors use (T/300)^beta
    std::vector<double> Jfit_coefs = {0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0};
    std::vector<double> Ffit_coefs = {0.0, 0.0, 0.0, 0.0};
    
    {
        // reaction 12:  CO2 + E => CO + O + E
        Janev_sum = 0.0;
        double k_f;
        Ffit_coefs = {
            28.3950215483559, -119.510009432235, 160.437496467960,
            -74.1425574458809};
        double Ffit_A = 1.50741518933862e-16;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2_ID] -= qdot;
    }

    {
        // reaction 13:  CO2 + E => CO2+ + E + E
        Janev_sum = 0.0;
        double k_f;
        Ffit_coefs = {
            65.5733184943082, -511.020275942352, 911.053057794656,
            -540.822599138513};
        double Ffit_A = 1.12781185638050e-15;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO2p_ID] += qdot;
        int rxntype = 2;
        double eexci = 13.8;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
        //amrex::Print()<<"CO2+ production:"<<wdot[CO2p_ID]<<"\t"<<13<<"\n";
    }

    {
        // reaction 18:  CO + E => C + O + E
        Janev_sum = 0.0;
        double k_f;
        Ffit_coefs = {
            16.3541660916144, -220.870249708484, 433.102467731706,
            -274.716891876346};
        double Ffit_A = 7.39606286753756e-15;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[CO_ID] -= qdot;
        wdot[C_ID] += qdot;
        int rxntype = 2;
        double eexci = 11.1;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 22:  AR + E => AR+ + E + E
        Janev_sum = 0.0;
        double k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -1951.40116682516, 9125.28938191463,  -19064.2751720459,
                22763.7019090326,  -16855.4371729404, 7904.47881474249,
                -2291.72587558317, 375.740475221758,  -26.6926001068518};
            double Jfit_A = 1.42900000000000e-14;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                -77.9510931269820, 533.060249822779, -2091.33857389100,
                1893.83145115372};
            double Ffit_A = 1.13390537981533e-12;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const double qf = k_f * (sc[E_ID] * sc[AR_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[AR_ID] -= qdot;
        wdot[ARp_ID] += qdot;
        int rxntype = 2;
        double eexci = 15.759;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 27:  H2 + E => 2 H + E
        Janev_sum = 0.0;
        double k_f;
        Ffit_coefs = {
            -3.28939221926166, -120.803981835609, 205.581141005113,
            -110.109964887129};
        double Ffit_A = 5.30049361173221e-15;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[H2_ID] -= qdot;
    }

    {
        // reaction 15:  H2 + E => 2 H + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -7.40026063650055,  1.84694302473460,   2.10477652143377,
            0.544578203292530,  -0.439269613111889, -0.190032372949602,
            0.0184361506891109, 0.0168253218120052, 0.00189898959859980};
        double Jfit_A = 7.88500000000000e-15;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f =
        std::min(Jfit_A * exp(Janev_sum) * 6.02214085774e23, 4.8e9);
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[H2_ID] -= qdot;
    }

    {
        // reaction 15:  H2 + E => 2 H + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -7.23752196537499,  1.84906543851286,   2.11123657619837,
            0.543622476041301,  -0.444326744928734, -0.191899221782850,
            0.0185268577418843, 0.0169571673407612, 0.00191449481252444};
        double Jfit_A = 6.69300000000000e-15;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f =
        std::min(Jfit_A * exp(Janev_sum) * 6.02214085774e23, 4.2e9);
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[H2_ID] -= qdot;
    }

    {
        // reaction 44:  H2 + E => E + H2+ + E
        Janev_sum = 0.0;
        double k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -1857.73884641047, 9635.17665706921,  -21676.8880597883,
                27257.3024281753,  -20922.3742900842, 10057.6678627611,
                -2964.28003383507, 490.871842318378,  -35.0368193361303};
            double Jfit_A = 5.28200000000000e-15;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                -35.6107865815828, 339.861096189597, -1989.12742271793,
                2664.29385213052};
            double Ffit_A = 3.21438910016961e-14;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2p_ID] += qdot;
        int rxntype = 2;
        double eexci = 15.40;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 45:  H2 + E => H + H-
        Janev_sum = 0.0;
        Jfit_coefs = {
            -5.24436622270143,  1.67023573577288,   1.39354323542527,
            0.244196472865026,  -0.312518441280522, -0.112078040187164,
            0.0161059139891788, 0.0110294867650971, 0.00118630502279277};
        double Jfit_A = 1.53600000000000e-18;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f =
        std::min(Jfit_A * exp(Janev_sum) * 6.02214085774e23, 9.0e5);
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H_ID] += qdot;
        wdot[Hn_ID] += qdot;
        wdot[H2_ID] -= qdot;
        int rxntype = 4;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 43:  E + O2 + M => O2- + M
        const double k_f = 1087985.3799976;
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[E_ID] * sc[O2_ID]);
        const double qr =
        Corr * k_f * exp(-(g_RT[0] + g_RT[9] - g_RT[10])) * (refC) * (0.0);
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[O2_ID] -= qdot;
        wdot[O2n_ID] += qdot;
        int rxntype = 4;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 44:  E + O + M => O- + M
        const double k_f = 36266.1793332534;
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[E_ID] * sc[O_ID]);
        const double qr =
        Corr * k_f * exp(-(g_RT[0] + g_RT[7] - g_RT[8])) * (refC) * (0.0);
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[O_ID] -= qdot;
        wdot[On_ID] += qdot;
        int rxntype = 4;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 50:  O2- + O2 => E + O2 + O2
        const double k_f = 1.31282668568;
        const double qf = k_f * (sc[O2_ID] * sc[O2n_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        int rxntype = 7;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 51:  O2- + H2O => E + O2 + H2O
        const double k_f = 3011070380;
        const double qf = k_f * (sc[O2n_ID] * sc[H2O_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        int rxntype = 7;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 52:  O2- + M => E + O2 + M
        const double k_f = 162597800.52 * exp((0.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[O2n_ID]);
        const double qr = Corr * k_f *
        exp(-(-g_RT[0] - g_RT[9] + g_RT[10])) *
        (refCinv) * (0.0);
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        int rxntype = 7;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 56:  H- + M => E + H + M
        const double k_f = 602214076;
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[Hn_ID]);
        const double qr = Corr * k_f *
        exp(-(-g_RT[0] - g_RT[1] + g_RT[2])) *
        (refCinv) * (0.0);
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[H_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        int rxntype = 7;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 84:  AR + ARe + M => AR2e + M
        const double k_f = 11967.8391799736;
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[38] * sc[39]);
        const double qr = Corr * k_f *
        exp(-(g_RT[38] + g_RT[39] - g_RT[42])) * (refC) *
        (0.0);
        const double qdot = qf - qr;
        wdot[AR_ID] -= qdot;
        wdot[ARe_ID] -= qdot;
        wdot[AR2e_ID] += qdot;
    }

    {
        // reaction 86:  AR2e + M => 2 AR + M
        const double k_f = 3011.07038;
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[42]);
        const double qr = Corr * k_f *
        exp(-(-2.000000 * g_RT[38] + g_RT[42])) *
        (refCinv) * (0.0);
        const double qdot = qf - qr;
        wdot[AR_ID] += 2.000000 * qdot;
        wdot[AR2e_ID] -= qdot;
    }

    {
        // reaction 87:  AR2e + M => AR + ARe + M
        const double k_f = 3011.07038;
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[42]);
        const double qr = Corr * k_f *
        exp(-(-g_RT[38] - g_RT[39] + g_RT[42])) *
        (refCinv) * (0.0);
        const double qdot = qf - qr;
        wdot[AR_ID] += qdot;
        wdot[ARe_ID] += qdot;
        wdot[AR2e_ID] -= qdot;
    }

    {
        // reaction 138:  H2v + M => H2 + M
        const double k_f = 60221.4076;
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[4]);
        const double qr = Corr * k_f * exp(-(-g_RT[3] + g_RT[4])) * (0.0);
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H2v_ID] -= qdot;
    }

    {
        // reaction 140:  AR+ + O- + M => AR + O + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[8] * sc[40]);
        const double qr =
        Corr * k_f * exp(-(-g_RT[7] + g_RT[8] - g_RT[38] + g_RT[40])) *
        (0.0);
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARp_ID] -= qdot;
    }

    {
        // reaction 143:  AR+ + O2- + M => AR + O2 + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[10] * sc[40]);
        const double qr =
        Corr * k_f * exp(-(-g_RT[9] + g_RT[10] - g_RT[38] + g_RT[40])) *
        (0.0);
        const double qdot = qf - qr;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARp_ID] -= qdot;
    }

    {
        // reaction 148:  AR+ + H- + M => AR + H + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[2] * sc[40]);
        const double qr =
        Corr * k_f * exp(-(-g_RT[1] + g_RT[2] - g_RT[38] + g_RT[40])) *
        (0.0);
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARp_ID] -= qdot;
    }

    {
        // reaction 151:  ARH+ + O- + M => AR + H + O + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[8] * sc[41]);
        const double qr =
        Corr * k_f *
        exp(-(-g_RT[1] - g_RT[7] + g_RT[8] - g_RT[38] + g_RT[41])) *
        (refCinv) * (0.0);
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARHp_ID] -= qdot;
    }

    {
        // reaction 154:  ARH+ + O2- + M => AR + H + O2 + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[10] * sc[41]);
        const double qr =
        Corr * k_f *
        exp(-(-g_RT[1] - g_RT[9] + g_RT[10] - g_RT[38] + g_RT[41])) *
        (refCinv) * (0.0);
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARHp_ID] -= qdot;
    }

    {
        // reaction 158:  H2+ + O- + M => H2O + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[5] * sc[8]);
        const double qr =
        Corr * k_f * exp(-(g_RT[5] + g_RT[8] - g_RT[24])) * (refC) * (0.0);
        const double qdot = qf - qr;
        wdot[H2p_ID] -= qdot;
        wdot[On_ID] -= qdot;
        wdot[H2O_ID] += qdot;
    }

    {
        // reaction 159:  H2+ + O- + M => H2 + O + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[5] * sc[8]);
        const double qr =
        Corr * k_f * exp(-(-g_RT[3] + g_RT[5] - g_RT[7] + g_RT[8])) * (0.0);
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
    }

    {
        // reaction 164:  O- + OH+ + M => O + OH + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[8] * sc[23]);
        const double qr =
        Corr * k_f * exp(-(-g_RT[7] + g_RT[8] - g_RT[22] + g_RT[23])) *
        (0.0);
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[OHp_ID] -= qdot;
    }

    {
        // reaction 167:  H2O+ + O- + M => H2O + O + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[8] * sc[25]);
        const double qr =
        Corr * k_f * exp(-(-g_RT[7] + g_RT[8] - g_RT[24] + g_RT[25])) *
        (0.0);
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
    }

    {
        // reaction 173:  H2+ + O2- + M => H2 + O2 + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[5] * sc[10]);
        const double qr = Corr * k_f *
        exp(-(-g_RT[3] + g_RT[5] - g_RT[9] + g_RT[10])) *
        (0.0);
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
    }

    {
        // reaction 179:  O2- + OH+ + M => O2 + OH + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[10] * sc[23]);
        const double qr =
        Corr * k_f * exp(-(-g_RT[9] + g_RT[10] - g_RT[22] + g_RT[23])) *
        (0.0);
        const double qdot = qf - qr;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[OHp_ID] -= qdot;
    }

    {
        // reaction 184:  H2O+ + O2- + M => H2O + O2 + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[10] * sc[25]);
        const double qr =
        Corr * k_f * exp(-(-g_RT[9] + g_RT[10] - g_RT[24] + g_RT[25])) *
        (0.0);
        const double qdot = qf - qr;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
    }

    {
        // reaction 192:  H- + H2+ + M => H + H2 + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[2] * sc[5]);
        const double qr =
        Corr * k_f * exp(-(-g_RT[1] + g_RT[2] - g_RT[3] + g_RT[5])) * (0.0);
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[H2p_ID] -= qdot;
    }

    {
        // reaction 195:  H- + OH+ + M => H + OH + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[2] * sc[23]);
        const double qr =
        Corr * k_f * exp(-(-g_RT[1] + g_RT[2] - g_RT[22] + g_RT[23])) *
        (0.0);
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[OHp_ID] -= qdot;
    }

    {
        // reaction 196:  H- + OH+ + M => H2O + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[2] * sc[23]);
        const double qr =
        Corr * k_f * exp(-(g_RT[2] + g_RT[23] - g_RT[24])) * (refC) * (0.0);
        const double qdot = qf - qr;
        wdot[Hn_ID] -= qdot;
        wdot[OHp_ID] -= qdot;
        wdot[H2O_ID] += qdot;
    }

    {
        // reaction 199:  H- + H2O+ + M => H + H2O + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[2] * sc[25]);
        const double qr =
        Corr * k_f * exp(-(-g_RT[1] + g_RT[2] - g_RT[24] + g_RT[25])) *
        (0.0);
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
    }

    {
        // reaction 202:  H- + H3O+ + M => H2 + H2O + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[2] * sc[26]);
        const double qr = Corr * k_f *
        exp(-(g_RT[2] - g_RT[3] - g_RT[24] + g_RT[26])) *
        (0.0);
        const double qdot = qf - qr;
        wdot[Hn_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[H2O_ID] += qdot;
        wdot[H3Op_ID] -= qdot;
    }

    {
        // reaction 203:  H- + H3O+ + M => H + H2 + OH + M
        const double k_f = 72532358666.5068 * exp((-2.5) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[2] * sc[26]);
        const double qr =
        Corr * k_f *
        exp(-(-g_RT[1] + g_RT[2] - g_RT[3] - g_RT[22] + g_RT[26])) *
        (refCinv) * (0.0);
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[OH_ID] += qdot;
        wdot[H3Op_ID] -= qdot;
    }

    {
        // reaction 261:  CO2 + CO2+ + M => C2O4+ + M
        const double k_f = 108798537.99976;
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[12] * sc[17]);
        const double qr = Corr * k_f *
        exp(-(g_RT[12] + g_RT[17] - g_RT[19])) * (refC) *
        (0.0);
        const double qdot = qf - qr;
        wdot[CO2_ID] -= qdot;
        wdot[CO2p_ID] -= qdot;
        wdot[C2O4p_ID] += qdot;
        //amrex::Print()<<"CO2+ production:"<<wdot[CO2p_ID]<<"\t"<<261<<"\t"<<sc[CO2_ID]<<"\t"<<sc[CO2p_ID]<<"\n";
    }

    {
        // reaction 263:  C2O4+ + CO + M => C2O3+ + CO2 + M
        const double k_f = 15231795319.9664;
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[11] * sc[19]);
        const double qr =
        Corr * k_f * exp(-(g_RT[11] - g_RT[12] - g_RT[18] + g_RT[19])) *
        (0.0);
        const double qdot = qf - qr;
        wdot[CO_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        wdot[C2O3p_ID] += qdot;
        wdot[C2O4p_ID] -= qdot;
    }

    {
        // reaction 264:  C2O4+ + M => CO2 + CO2+ + M
        const double k_f = 6022.14076;
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[19]);
        const double qr = Corr * k_f *
        exp(-(-g_RT[12] - g_RT[17] + g_RT[19])) *
        (refCinv) * (0.0);
        const double qdot = qf - qr;
        wdot[CO2_ID] += qdot;
        wdot[CO2p_ID] += qdot;
        wdot[C2O4p_ID] -= qdot;
        //amrex::Print()<<"CO2+ production:"<<wdot[CO2p_ID]<<"\t"<<264<<"\n";
    }

    {
        // reaction 270:  CO2 + O- + M => CO3- + M
        const double k_f = 32639561.399928;
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[8] * sc[12]);
        const double qr =
        Corr * k_f * exp(-(g_RT[8] + g_RT[12] - g_RT[20])) * (refC) * (0.0);
        const double qdot = qf - qr;
        wdot[On_ID] -= qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO3n_ID] += qdot;
    }

    {
        // reaction 274:  2 O + AR => O2 + AR
        const double k_f =
        18.894679432625 * exp(-(-899.751398458839) * invT);
        const double qf = k_f * ((sc[7] * sc[7]) * sc[38]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= 2.000000 * qdot;
        wdot[O2_ID] += qdot;
    }

    {
        // reaction 275:  2 O + O2 => O2 + O2
        const double k_f =
        18.894679432625 * exp(-(-899.751398458839) * invT);
        const double qf = k_f * ((sc[7] * sc[7]) * sc[9]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= 2.000000 * qdot;
        wdot[O2_ID] += qdot;
    }

    {
        // reaction 276:  2 O + H2 => O2 + H2
        const double k_f =
        47.1460331332294 * exp(-(-899.751398458839) * invT);
        const double qf = k_f * (sc[3] * (sc[7] * sc[7]));
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= 2.000000 * qdot;
        wdot[O2_ID] += qdot;
    }

    {
        // reaction 277:  2 O + H2O => O2 + H2O
        const double k_f =
        94.2920662664588 * exp(-(-899.751398458839) * invT);
        const double qf = k_f * ((sc[7] * sc[7]) * sc[24]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= 2.000000 * qdot;
        wdot[O2_ID] += qdot;
    }

    {
        // reaction 278:  H + O + AR => OH + AR
        const double k_f = 5439.92689998801 * exp((-1) * logT);
        const double qf = k_f * (sc[1] * sc[7] * sc[38]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[O_ID] -= qdot;
        wdot[OH_ID] += qdot;
    }

    {
        // reaction 279:  H + O + O2 => OH + O2
        const double k_f = 5439.92689998801 * exp((-1) * logT);
        const double qf = k_f * (sc[1] * sc[7] * sc[9]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[O_ID] -= qdot;
        wdot[OH_ID] += qdot;
    }

    {
        // reaction 280:  H + O + H2 => OH + H2
        const double k_f = 10879.853799976 * exp((-1) * logT);
        const double qf = k_f * (sc[1] * sc[3] * sc[7]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[O_ID] -= qdot;
        wdot[OH_ID] += qdot;
    }

    {
        // reaction 281:  H + O + H2O => OH + H2O
        const double k_f = 27199.63449994 * exp((-1) * logT);
        const double qf = k_f * (sc[1] * sc[7] * sc[24]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[O_ID] -= qdot;
        wdot[OH_ID] += qdot;
    }

    {
        // reaction 283:  O2 + H => 2 O + H
        const double k_f = 3613284456 * exp(-(52294.2759104265) * invT);
        const double qf = k_f * (sc[1] * sc[9]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += 2.000000 * qdot;
        wdot[O2_ID] -= qdot;
    }

    {
        // reaction 284:  2 H + AR => H2 + AR
        const double k_f = 7253.23586665068 * exp((-1) * logT);
        const double qf = k_f * ((sc[1] * sc[1]) * sc[38]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= 2.000000 * qdot;
        wdot[H2_ID] += qdot;
    }

    {
        // reaction 285:  2 H + O2 => H2 + O2
        const double k_f = 7253.23586665068 * exp((-1) * logT);
        const double qf = k_f * ((sc[1] * sc[1]) * sc[9]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= 2.000000 * qdot;
        wdot[H2_ID] += qdot;
    }

    {
        // reaction 286:  2 H + H2 => H2 + H2
        const double k_f = 14506.4717333014 * exp((-1) * logT);
        const double qf = k_f * ((sc[1] * sc[1]) * sc[3]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= 2.000000 * qdot;
        wdot[H2_ID] += qdot;
    }

    {
        // reaction 287:  2 H + H2O => H2 + H2O
        const double k_f = 33364.8849865931 * exp((-1) * logT);
        const double qf = k_f * ((sc[1] * sc[1]) * sc[24]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= 2.000000 * qdot;
        wdot[H2_ID] += qdot;
    }

    {
        // reaction 288:  H2 + H => 2 H + H
        const double k_f =
        281233973492 * exp((-1) * logT - (54994.0333224688) * invT);
        const double qf = k_f * (sc[1] * sc[3]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[H2_ID] -= qdot;
    }

    {
        // reaction 289:  OH + H => H + O + H
        const double k_f = 3613284456 * exp(-(50894.3271461577) * invT);
        const double qf = k_f * (sc[1] * sc[22]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += qdot;
        wdot[OH_ID] -= qdot;
    }

    {
        // reaction 291:  H + OH + AR => H2O + AR
        const double k_f = 290129.434666027 * exp((-2.6) * logT);
        const double qf = k_f * (sc[1] * sc[22] * sc[38]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[OH_ID] -= qdot;
        wdot[H2O_ID] += qdot;
    }

    {
        // reaction 292:  H + OH + O2 => H2O + O2
        const double k_f = 290129.434666027 * exp((-2.6) * logT);
        const double qf = k_f * (sc[1] * sc[9] * sc[22]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[OH_ID] -= qdot;
        wdot[H2O_ID] += qdot;
    }

    {
        // reaction 293:  H + OH + H2 => H2O + H2
        const double k_f = 652791.227998561 * exp((-2.6) * logT);
        const double qf = k_f * (sc[1] * sc[3] * sc[22]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[OH_ID] -= qdot;
        wdot[H2O_ID] += qdot;
    }

    {
        // reaction 294:  H + OH + H2O => H2O + H2O
        const double k_f = 1447020.55539681 * exp((-2.6) * logT);
        const double qf = k_f * (sc[1] * sc[22] * sc[24]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[OH_ID] -= qdot;
        wdot[H2O_ID] += qdot;
    }

    {
        // reaction 295:  H2O + H => H + OH + H
        const double k_f = 3492841640.8 * exp(-(52894.1101760657) * invT);
        const double qf = k_f * (sc[1] * sc[24]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[OH_ID] += qdot;
        wdot[H2O_ID] -= qdot;
    }

    {
        // reaction 297:  O2 + H2 => 2 O + H2
        const double k_f = 3613284456 * exp(-(52294.2759104265) * invT);
        const double qf = k_f * (sc[3] * sc[9]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += 2.000000 * qdot;
        wdot[O2_ID] -= qdot;
    }

    {
        // reaction 300:  OH + H2 => H + O + H2
        const double k_f = 3613284456 * exp(-(50894.3271461577) * invT);
        const double qf = k_f * (sc[3] * sc[22]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += qdot;
        wdot[OH_ID] -= qdot;
    }

    {
        // reaction 301:  H2O + H2 => H + OH + H2
        const double k_f = 3492841640.8 * exp(-(52894.1101760657) * invT);
        const double qf = k_f * (sc[3] * sc[24]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[OH_ID] += qdot;
        wdot[H2O_ID] -= qdot;
    }

    {
        // reaction 308:  CH3 + H + M => CH4 + M
        const double k_f = 109161199.793093 * exp((-1.8) * logT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[1] * sc[31]);
        const double qr =
        Corr * k_f * exp(-(g_RT[1] + g_RT[31] - g_RT[33])) * (refC) * (0.0);
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[CH3_ID] -= qdot;
        wdot[CH4_ID] += qdot;
    }

    {
        // reaction 324:  CH3 + OH + M => CH3OH + M
        const double k_f =
        1338222017.39705 * exp(-(1279.6799811414) * invT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[22] * sc[31]);
        const double qr = Corr * k_f *
        exp(-(g_RT[22] + g_RT[31] - g_RT[36])) * (refC) *
        (0.0);
        const double qdot = qf - qr;
        wdot[OH_ID] -= qdot;
        wdot[CH3_ID] -= qdot;
        wdot[CH3OH_ID] += qdot;
    }

    {
        // reaction 334:  CO + H + M => HCO + M
        const double k_f =
        721.696968731742 * exp(-(841.881481891296) * invT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[1] * sc[11]);
        const double qr =
        Corr * k_f * exp(-(g_RT[1] + g_RT[11] - g_RT[27])) * (refC) * (0.0);
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[CO_ID] -= qdot;
        wdot[HCO_ID] += qdot;
    }

    {
        // reaction 374:  CO2 + M => CO + O + M
        const double k_f = 235465703.716 * exp(-(49424.4312653422) * invT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[12]);
        const double qr = Corr * k_f *
        exp(-(-g_RT[7] - g_RT[11] + g_RT[12])) *
        (refCinv) * (0.0);
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2_ID] -= qdot;
    }

    {
        // reaction 376:  CO + O + M => CO2 + M
        const double k_f =
        297.382670532678 * exp(-(1509.64999741416) * invT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[7] * sc[11]);
        const double qr =
        Corr * k_f * exp(-(g_RT[7] + g_RT[11] - g_RT[12])) * (refC) * (0.0);
        const double qdot = qf - qr;
        wdot[O_ID] -= qdot;
        wdot[CO_ID] -= qdot;
        wdot[CO2_ID] += qdot;
    }

    {
        // reaction 378:  C + O + M => CO + M
        const double k_f =
        7760962.37731622 * exp((-3.08) * logT - (2114.01321304563) * invT);
        const double Corr = mixture;
        const double qf = Corr * k_f * (sc[7] * sc[21]);
        const double qr =
        Corr * k_f * exp(-(g_RT[7] - g_RT[11] + g_RT[21])) * (refC) * (0.0);
        const double qdot = qf - qr;
        wdot[O_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[C_ID] -= qdot;
    }

    {
        // reaction 0:  CO2 + E => CO2 + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.74567823745404,   0.489892009971784,    0.475497844841483,
            -0.0828111482837000, -0.0869447627283896,  0.00364885418266187,
            0.00738742037568563, 0.000815768871163232, -1.80446057296271e-05};
        double Jfit_A = 3.07100000000000e-13;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f =
        std::min(Jfit_A * exp(Janev_sum) * 6.02214085774e23, 1.0e12);
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        int rxntype = 5;
        double eexci = 0.0;
        int elidx = CO2_ID;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: excitational losses for CO2(e1)
        // reaction 0:  CO2 + E => CO2 + E
        Janev_sum = 0.0;
        double k_f;
        Ffit_coefs = {
            13.9944258667849, -44.2367004233860, 35.3755884340306,
            -12.2196503543953};
        double Ffit_A = 4.66040157646352e-16;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        int rxntype = 2;
        double eexci = 7.0;
        int elidx = CO2_ID;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: excitational losses for CO2(e2)
        // reaction 0:  CO2 + E => CO2 + E
        Janev_sum = 0.0;
        double k_f;
        Ffit_coefs = {
            -12.6852242013963, -39.8289588750467, 75.9272218485778,
            -44.8968368578706};
        double Ffit_A = 1.41992889938866e-13;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        int rxntype = 2;
        double eexci = 10.5;
        int elidx = CO2_ID;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 1:  CO2 + E => CO2v1 + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -0.502839945828532, 0.694016395482428,    0.190052560381054,
            -0.336722426572758, -0.112995969891858,   0.0625050243076499,
            0.0180812366062628, -0.00406798152564209, -0.00117479961936418};
        double Jfit_A = 9.13200000000000e-15;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO2v1_ID] += qdot;
        int rxntype = 1;
        double eexci = 0.083;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: CO2(020)
        // reaction 2:  CO2 + E => CO2v2 + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.08379968039393,  1.78557074369457,    0.0320873806354341,
            -0.785497270396418, -0.0745881449840463, 0.173672656884839,
            0.0242491592490644, -0.0137085592778002, -0.00283607126964260};
        double Jfit_A = 4.31900000000000e-15;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO2v2_ID] += qdot;
        int rxntype = 1;
        double eexci = 0.17;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: CO2(100)
        // reaction 3:  CO2 + E => CO2v2 + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.74710336851205,   5.28866024293583,    -3.06722524904491,
            -1.26813972992265,   0.765286490697251,   0.220019094021629,
            -0.0660599223620025, -0.0174513395327540, -0.000100232840708316};
        double Jfit_A = 2.10800000000000e-15;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO2v2_ID] += qdot;
        int rxntype = 1;
        double eexci = 0.17;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 4:  CO2 + E => CO2v3 + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -0.0161587259013659, 0.192761580775350,    -0.233095798484862,
            -0.0504138326883387, 0.0729053506151256,   0.0268158476130099,
            -0.0130233553843901, -0.00234706991825410, 0.000208569444277222};
        double Jfit_A = 5.87800000000000e-15;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO2v3_ID] += qdot;
        int rxntype = 1;
        double eexci = 0.291;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 5:  CO2 + E => CO2v4 + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.26399041657649,  4.01593315343776,  -3.86564987870520,
            0.630843859313518,  0.918455783525699, -0.466408614313122,
            -0.122041546730277, 0.133320282875472, -0.0249335094184505};
        double Jfit_A = 1.22800000000000e-15;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO2v4_ID] += qdot;
        int rxntype = 1;
        double eexci = 0.252;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 6:  CO2 + E => CO2v4 + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -0.430350115379283,  1.88554119808260,   -2.10311501560593,
            0.205835857388545,   0.349603651488484,  -0.113940765190612,
            -0.0356282206715671, 0.0288224520101671, -0.00505168357993776};
        double Jfit_A = 6.65200000000000e-16;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO2v4_ID] += qdot;
        int rxntype = 1;
        double eexci = 0.339;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 7:  CO2 + E => CO2v4 + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.26866503084565,  4.16951503938324,  -4.27011214236418,
            0.528779192020846,  1.67719403845492,  -0.767308824575880,
            -0.423310590818580, 0.369354928743402, -0.0691523136873366};
        double Jfit_A = 1.02300000000000e-15;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO2v4_ID] += qdot;
        int rxntype = 1;
        double eexci = 0.339;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 8:  CO2 + E => CO2v4 + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.16752755372413,  3.97264512328113,  -4.02945519654257,
            0.673891717708656,  1.00483069462086,  -0.508781573502289,
            -0.133736912597309, 0.144832888470466, -0.0269944163072462};
        double Jfit_A = 3.90400000000000e-16;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO2v4_ID] += qdot;
        int rxntype = 1;
        double eexci = 0.422;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 9:  CO2 + E => CO2v4 + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.57218058544877,  4.53507634411035,  -3.90997647678875,
            0.335368207124721,  1.05569484438475,  -0.410009273039837,
            -0.159898555151304, 0.129957825662046, -0.0219563391933853};
        double Jfit_A = 2.89300000000000e-16;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO2v4_ID] += qdot;
        int rxntype = 1;
        double eexci = 0.5;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 10:  CO2 + E => CO2v4 + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -0.803462147557306,  3.31431359636670,   -3.80327983310172,
            0.791648348030451,   0.692211498451406,  -0.342472582689845,
            -0.0688634709416767, 0.0708789327184151, -0.0117253943723200};
        double Jfit_A = 1.68100000000000e-16;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO2v4_ID] += qdot;
        int rxntype = 1;
        double eexci = 0.505;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 11:  CO2 + E => CO2v4 + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.13341741015961,  3.84235032466761,  -3.89843865763089,
            0.695260507371384,  0.918161189362339, -0.486298345254785,
            -0.117448184645462, 0.135925105891503, -0.0258281763494971};
        double Jfit_A = 5.46700000000000e-16;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CO2v4_ID] += qdot;
        int rxntype = 1;
        double eexci = 2.5;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 16:  CO2 + E => CO + O-
        Janev_sum = 0.0;
        double k_f;
        Ffit_coefs = {
            5.68537810168761, -15.1912011431635, 6.86995641002507,
            -1.21745244042248};
        double Ffit_A = 4.23133262631561e-18;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[CO2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[On_ID] += qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        int rxntype = 4;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 15:  CO + E => CO + E
        const double k_f = 0;
        const double qf = k_f * (sc[0] * sc[11]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[CO_ID] -= qdot;
        wdot[CO_ID] += qdot;
    }

    {
        // reaction 20:  AR + E => AR + E
        Janev_sum = 0.0;
        Jfit_coefs = {
            -2.62130096848181,    1.55888667906962,     0.0350711485300335,
            0.00439515704875221,  -0.00138205864597430, -0.0293026895863627,
            -0.00284458054782676, 0.00241786232297547,  0.000390470745241025};
        double Jfit_A = 2.33300000000000e-13;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = std::min(
            Jfit_A * exp(Janev_sum) * 6.02214085774e23,
            3.0e-13 * 6.02214085774e23);
        const double qf = k_f * (sc[E_ID] * sc[AR_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[AR_ID] -= qdot;
        wdot[AR_ID] += qdot;
        int rxntype = 5;
        double eexci = 0.0;
        int elidx = AR_ID;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 21:  AR + E => ARe + E
        Janev_sum = 0.0;
        double k_f;
        Ffit_coefs = {
            -17.8244315740958, -55.1739383730274, 102.978958246859,
            -55.9966048149755};
        double Ffit_A = 9.18800465289008e-14;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[AR_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[AR_ID] -= qdot;
        wdot[ARe_ID] += qdot;
        int rxntype = 1;
        double eexci = 11.55;
        int elidx = AR_ID;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 23:  E + H2 => E + H2
        Janev_sum = 0.0;
        Jfit_coefs = {
            -0.158926451167132,    0.393305952760301,     -0.214034386018838,
            -0.0313291759745509,   0.0107225576024030,    0.00170840526110980,
            -0.000413389366430107, -5.14813979990267e-05, 3.19815688057637e-06};
        double Jfit_A = 1.27600000000000e-13;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 5;
        double eexci = 0.0;
        int elidx = H2_ID;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 24:  E + H2 => E + H2v
        Janev_sum = 0.0;
        Jfit_coefs = {
            -0.804296492019153,  1.53571652805659,     -0.697157741030016,
            -0.0743282818188661, 0.0184933104510412,   0.0226678142666104,
            0.00107857037424512, -0.00120726454339529, -0.000496070060173857};
        double Jfit_A = 3.99600000000000e-15;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2v_ID] += qdot;
        int rxntype = 1;
        double eexci = 0.5;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy loss to H2(v2)
        // reaction 25:  E + H2 => E + H2v
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.29886441235163,  2.36736721872670,  -1.56670710235572,
            0.270073946043366,  0.488592533709920, -0.316604030640732,
            -0.103629264418870, 0.117047458324374, -0.0224406222651364};
        double Jfit_A = 3.12000000000000e-16;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2v_ID] += qdot;
        int rxntype = 1;
        double eexci = 1.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy loss to H2(v3)
        // reaction 26:  E + H2 => E + H2v
        Janev_sum = 0.0;
        Jfit_coefs = {
            -1.32007485064830,   2.57319732174652,   -1.75221503040112,
            0.357811634017136,   0.334804775754479,  -0.248830468289447,
            -0.0422357401405263, 0.0656821644773556, -0.0125518656838086};
        double Jfit_A = 3.00800000000000e-17;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2v_ID] += qdot;
        int rxntype = 1;
        double eexci = 1.5;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H2(b3Su) reaction 28:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        Ffit_coefs = {
            -2.75446684304778, -50.2291122564490, 73.9927260539435,
            -36.8788754298209};
        double Ffit_A = 9.72983861974428e-15;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 7.93;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H2(B1Su) reaction 29:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        Ffit_coefs = {
            9.73715263713406, -169.759074219013, 302.650832488494,
            -173.548526278283};
        double Ffit_A = 7.08804984129738e-15;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 11.40;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H2(a3Sg) reaction 30:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -1979.60623192343, 10192.8890713498,  -22690.2521657955,
                28370.3653882692,  -21785.5023661081, 10534.0833206771,
                -3137.20812255726, 526.924450354538,  -38.2654901887211};
            double Jfit_A = 6.09300000000000e-16;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                -46.8536837163814, 444.211477664037, -1919.21546688701,
                2414.84708251725};
            double Ffit_A = 4.14543774825509e-15;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 11.72;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H2(c3Pu) reaction 31:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -1296.16503964429, 6693.57516651356,  -14951.3300656748,
                18738.8841595928,  -14406.8434523910, 6967.71435350915,
                -2073.99857605575, 347.983786681218,  -25.2354823022991};
            double Jfit_A = 1.36100000000000e-15;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                -58.4021251309784, 537.327503682081, -2190.50724519666,
                2749.52795344136};
            double Ffit_A = 1.49575563289077e-14;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 11.72;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H2(C1Pu) reaction 32:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -1722.95883904265, 8840.57300247299,  -19662.3928799025,
                24571.8767044338,  -18862.3663163976, 9120.30925661930,
                -2717.10515978014, 456.709922480483,  -33.2048491354566};
            double Jfit_A = 3.69800000000000e-15;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                -45.3587148056958, 333.894802125744, -1528.57489939913,
                1996.74422025009};
            double Ffit_A = 4.58123944995702e-14;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 12.40;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H2(EF1Sg) reaction 33:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -771.411775136490, 3875.96373773952,  -8439.31791893315,
                10225.2380718685,  -7518.89528497335, 3439.65743931899,
                -957.637795763425, 148.553978141686,  -9.83730600520898};
            double Jfit_A = 1.29200000000000e-15;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                -1.03257453231928, 72.4645068119013, -768.779715373533,
                1149.53632388452};
            double Ffit_A = 1.26310016561220e-15;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 12.40;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H2(e3Su) reaction 34:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -1150.55546835788, 6037.37189716131,  -13727.2170157118,
                17465.5451487199,  -13579.9508418143, 6617.34753254167,
                -1977.86041355528, 332.238089036930,  -24.0593794249128};
            double Jfit_A = 1.88600000000000e-16;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                -27.0484052355503, 300.601693287326, -1430.69852408452,
                1763.93654818736};
            double Ffit_A = 4.76340471972864e-16;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 13.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H2(B'1Su) reaction 35:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -2479.77794589974, 12868.7484710841,  -28939.2708170773,
                36569.2236606705,  -28384.4187726694, 13874.3691698582,
                -4177.29429555056, 709.329977389421,  -52.0768851998328};
            double Jfit_A = 4.33900000000000e-16;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                -77.3422255188623, 592.975423063700, -2476.31969986974,
                3079.38992638324};
            double Ffit_A = 2.30804716936867e-14;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 13.8;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H2(D1Pu) reaction 36:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -1439.85072838354, 7501.07414409844,  -16959.9526218994,
                21451.2488398011,  -16577.7782773135, 8030.15822703666,
                -2386.53547420628, 398.746313507096,  -28.7308290010795};
            double Jfit_A = 3.03500000000000e-16;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                -32.0167399284743, 259.526573297002, -1445.35098524570,
                1911.87694646186};
            double Ffit_A = 1.85948818153970e-15;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 14.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H2(B"1Su) reaction 37:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -1511.68343766125, 7853.26058824394,  -17710.3365190378,
                22336.3223503825,  -17205.4286954925, 8303.74570974656,
                -2458.02279114341, 408.947491419507,  -29.3345576577653};
            double Jfit_A = 4.31800000000000e-18;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                -30.1156005028426, 269.029182721723, -1564.93411061652,
                2088.48574093658};
            double Ffit_A = 2.18081518255907e-17;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 14.6;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H2(D'1Pu) reaction 38:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -1623.83862891674, 8436.33548817704,  -19010.5486957234,
                23949.7480501913,  -18423.4444914014, 8877.64829120409,
                -2623.24137908025, 435.579097522363,  -31.1783913722693};
            double Jfit_A = 7.01500000000000e-17;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                -33.4626924276693, 315.067335486014, -1763.06821892842,
                2328.44773677504};
            double Ffit_A = 3.69314852331613e-16;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 14.6;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H(1S) + H(2P) reaction 39:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -1040.41862142214, 5184.68998652417,  -11244.7072600707,
                13561.7679057104,  -9911.01285847446, 4499.27707708176,
                -1241.34744203749, 190.568117411283,  -12.4697342945185};
            double Jfit_A = 5.59500000000000e-16;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                5.39236346695952, 24.7748709446991, -971.786226819569,
                1598.10574023187};
            double Ffit_A = 5.22063789288094e-16;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 14.68;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H(1S) + H(2S) reaction 40:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        if (TeeV < 8.0)
        {
            Jfit_coefs = {
                -1869.33383446626, 9791.59892240395,  -22211.3783084400,
                28178.4213835468,  -21840.9619443652, 10609.0419753568,
                -3160.96346234453, 529.332662415895,  -38.2157977421070};
            double Jfit_A = 3.91200000000000e-16;
            for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
            k_f = Jfit_A * exp(Janev_sum) * 6.02214085774e23;
        } else
        {
            Ffit_coefs = {
                -56.1489806582272, 498.226575245362, -2285.71522518037,
                2826.91360708509};
            double Ffit_A = 5.17082365264867e-15;
            for (int j = 0; j < 4; j++)
                Janev_sum += Ffit_coefs[j] * invTe_pow[j];
            k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        }
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 14.68;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H(1S) + H(3) reaction 41:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        Ffit_coefs = {
            -10.4929158785148, -166.837354469033, 1287.16071993509,
            -4035.63968276149};
        double Ffit_A = 6.57526974996034e-16;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 16.57;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H(1S) + H(4) reaction 42:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        Ffit_coefs = {
            23.7861184571555, -579.496649338752, 3175.06620048815,
            -7295.94567126596};
        double Ffit_A = 5.39806251555896e-17;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 17.22;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: Accounting for electron energy losses only:  E + H2
        // => E + H(1S) + H(5) reaction 43:  E + H2 => E + H2
        Janev_sum = 0.0;
        double k_f;
        Ffit_coefs = {
            16.1202627782974, -533.444622444002, 3066.28946754340,
            -7084.75235705944};
        double Ffit_A = 2.18762463103676e-17;
        for (int j = 0; j < 4; j++) Janev_sum += Ffit_coefs[j] * invTe_pow[j];
        k_f = Ffit_A * exp(Janev_sum) * 6.02214085774e23;
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 1;
        double eexci = 17.53;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 45:  H2 + E => H + H-
        Janev_sum = 0.0;
        Jfit_coefs = {
            -5.24436622270143,  1.67023573577288,   1.39354323542527,
            0.244196472865026,  -0.312518441280522, -0.112078040187164,
            0.0161059139891788, 0.0110294867650971, 0.00118630502279277};
        double Jfit_A = 1.53600000000000e-18;
        for (int j = 0; j < 9; j++) Janev_sum += Jfit_coefs[j] * Te_pow[j];
        const double k_f =
        std::min(Jfit_A * exp(Janev_sum) * 6.02214085774e23, 9.0e5);
        const double qf = k_f * (sc[E_ID] * sc[H2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H_ID] += qdot;
        wdot[Hn_ID] += qdot;
        wdot[H2_ID] -= qdot;
        int rxntype = 4;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 45:  O + O- => E + O2
        const double k_f = 90332111.4;
        const double qf = k_f * (sc[O_ID] * sc[On_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[O_ID] -= qdot;
        wdot[On_ID] -= qdot;
        wdot[O2_ID] += qdot;
        int rxntype = 7;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 46:  H + O- => E + OH
        const double k_f = 301107038;
        const double qf = k_f * (sc[H_ID] * sc[On_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[H_ID] -= qdot;
        wdot[On_ID] -= qdot;
        wdot[OH_ID] += qdot;
        int rxntype = 7;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 47:  H2 + O- => E + H2O
        const double k_f = 404687859.072;
        const double qf = k_f * (sc[H2_ID] * sc[On_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[On_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        int rxntype = 7;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 48:  C + O- => CO + E
        const double k_f = 301107038;
        const double qf = k_f * (sc[On_ID] * sc[C_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[On_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[C_ID] -= qdot;
        int rxntype = 7;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 49:  CO + O- => CO2 + E
        const double k_f = 391439149.4;
        const double qf = k_f * (sc[On_ID] * sc[CO_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[On_ID] -= qdot;
        wdot[CO_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        int rxntype = 7;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 53:  H + H- => E + H2
        const double k_f = 782878298.8;
        const double qf = k_f * (sc[H_ID] * sc[Hn_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[H_ID] -= qdot;
        wdot[Hn_ID] -= qdot;
        wdot[H2_ID] += qdot;
        int rxntype = 7;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 54:  H- + O => E + OH
        const double k_f = 602214076;
        const double qf = k_f * (sc[Hn_ID] * sc[O_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[O_ID] -= qdot;
        wdot[OH_ID] += qdot;
        int rxntype = 7;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 55:  H- + OH => E + H2O
        const double k_f = 60221407.6;
        const double qf = k_f * (sc[Hn_ID] * sc[OH_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[OH_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        int rxntype = 7;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 57:  CH3 + H- => CH4 + E
        const double k_f = 602214076;
        const double qf = k_f * (sc[Hn_ID] * sc[CH3_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[CH3_ID] -= qdot;
        wdot[CH4_ID] += qdot;
        int rxntype = 7;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 58:  CH2 + H- => CH3 + E
        const double k_f = 602214076;
        const double qf = k_f * (sc[Hn_ID] * sc[CH2_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[CH2_ID] -= qdot;
        wdot[CH3_ID] += qdot;
        int rxntype = 7;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 59:  CO + CO3- => 2 CO2 + E
        const double k_f = 33.12177418;
        const double qf = k_f * (sc[CO_ID] * sc[CO3n_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[CO_ID] -= qdot;
        wdot[CO2_ID] += 2.000000 * qdot;
        wdot[CO3n_ID] -= qdot;
        int rxntype = 7;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 60:  CO2+ + E => CO + O
        const double k_f = 1294945000000000 * pow(Te, -0.5) / tc[1];
        const double qf = k_f * (sc[E_ID] * sc[CO2p_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[O_ID] += qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2p_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
        //amrex::Print()<<"CO2+ production:"<<wdot[CO2p_ID]<<"\t"<<60<<"\n";
    }

    {
        // reaction 61:  CO2+ + E => C + O2
        const double k_f = 9998180000000 * pow(Te, -0.4);
        const double qf = k_f * (sc[E_ID] * sc[CO2p_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[O2_ID] += qdot;
        wdot[CO2p_ID] -= qdot;
        wdot[C_ID] += qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
        //amrex::Print()<<"CO2+ production:"<<wdot[CO2p_ID]<<"\t"<<61<<"\n";
    }

    {
        // reaction 62:  C2O3+ + E => CO + CO2
        const double k_f = 22766940000000 * pow(Te, -0.7);
        const double qf = k_f * (sc[E_ID] * sc[C2O3p_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2_ID] += qdot;
        wdot[C2O3p_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 63:  C2O4+ + E => 2 CO2
        const double k_f = 1294945000000000 * pow(Te, -0.5) / tc[1];
        const double qf = k_f * (sc[E_ID] * sc[C2O4p_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[CO2_ID] += 2.000000 * qdot;
        wdot[C2O4p_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: estimate based on E + N2+ => 2N
        // reaction 91:  E + H2+ => 2 H
        const double k_f = 54207000000 * pow(Te / 300.0, -0.39);
        const double qf = k_f * (sc[E_ID] * sc[H2p_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[H2p_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: estimate based on E + N2+ => 2N
        // reaction 93:  E + H3+ => 3 H
        const double k_f = 54207000000 * pow(Te / 300.0, -0.39);
        const double qf = k_f * (sc[E_ID] * sc[H3p_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H_ID] += 3.000000 * qdot;
        wdot[H3p_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 95:  CH3+ + E => CH2 + H
        const double k_f = 13551750000 * pow(Te / 300, -0.5);
        const double qf = k_f * (sc[E_ID] * sc[CH3p_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H_ID] += qdot;
        wdot[CH2_ID] += qdot;
        wdot[CH3p_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 96:  CH3+ + E => C + H + H2
        const double k_f = 10178870000 * pow(Te / 300, -0.5);
        const double qf = k_f * (sc[E_ID] * sc[CH3p_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H_ID] += qdot;
        wdot[H2_ID] += qdot;
        wdot[C_ID] += qdot;
        wdot[CH3p_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 69:  CH2+ + E => C + H2
        const double k_f = 2903086000 * pow(Te / 300, -0.5);
        const double qf = k_f * (sc[E_ID] * sc[CH2p_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[C_ID] += qdot;
        wdot[CH2p_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 70:  CH2+ + E => C + 2 H
        const double k_f = 15238190000 * pow(Te / 300, -0.5);
        const double qf = k_f * (sc[E_ID] * sc[CH2p_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[C_ID] += qdot;
        wdot[CH2p_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 71:  CH+ + E => C + H
        const double k_f = 19454290000 * pow(Te / 300, -0.42);
        const double qf = k_f * (sc[E_ID] * sc[CHp_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H_ID] += qdot;
        wdot[C_ID] += qdot;
        wdot[CHp_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 72:  E + H3O+ => H + H2O
        const double k_f = 15057500000 * pow(TeeV, -0.7);
        const double qf = k_f * (sc[E_ID] * sc[H3Op_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H_ID] += qdot;
        wdot[H2O_ID] += qdot;
        wdot[H3Op_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 73:  E + H3O+ => H2 + OH
        const double k_f = 8432200000 * pow(TeeV, -0.7);
        const double qf = k_f * (sc[E_ID] * sc[H3Op_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[OH_ID] += qdot;
        wdot[H3Op_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 74:  E + H3O+ => 2 H + OH
        const double k_f = 36138000000 * pow(TeeV, -0.7);
        const double qf = k_f * (sc[E_ID] * sc[H3Op_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[OH_ID] += qdot;
        wdot[H3Op_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 75:  E + H3O+ => H + H2 + O
        const double k_f = 782990000 * pow(TeeV, -0.7);
        const double qf = k_f * (sc[E_ID] * sc[H3Op_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H_ID] += qdot;
        wdot[H2_ID] += qdot;
        wdot[O_ID] += qdot;
        wdot[H3Op_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 76:  E + H2O+ => H + OH
        const double k_f = 897427000000 * pow(Te, -0.5);
        const double qf = k_f * (sc[E_ID] * sc[H2Op_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H_ID] += qdot;
        wdot[OH_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 77:  E + H2O+ => H2 + O
        const double k_f = 286092500000 * pow(Te, -0.5);
        const double qf = k_f * (sc[E_ID] * sc[H2Op_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[O_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 78:  E + H2O+ => 2 H + O
        const double k_f = 3180144000000 * pow(Te, -0.5);
        const double qf = k_f * (sc[E_ID] * sc[H2Op_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[O_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 79:  E + OH+ => H + O
        const double k_f = 391495000000 * pow(Te, -0.5);
        const double qf = k_f * (sc[E_ID] * sc[OHp_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H_ID] += qdot;
        wdot[O_ID] += qdot;
        wdot[OHp_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction comment: estimate based on E + N2+ => 2N
        // reaction 111:  ARH+ + E => AR + H
        const double k_f = 54207000000 * pow(Te / 300.0, -0.39);
        const double qf = k_f * (sc[E_ID] * sc[ARHp_ID]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] -= qdot;
        wdot[H_ID] += qdot;
        wdot[AR_ID] += qdot;
        wdot[ARHp_ID] -= qdot;
        int rxntype = 3;
        double eexci = 0.0;
        int elidx = 0;
        comp_ener_exch(
            qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
    }

    {
        // reaction 81:  AR2e + ARe => 2 AR + AR+ + E
        const double k_f = 301107038;
        const double qf = k_f * (sc[39] * sc[42]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[AR_ID] += 2.000000 * qdot;
        wdot[ARe_ID] -= qdot;
        wdot[ARp_ID] += qdot;
        wdot[AR2e_ID] -= qdot;
    }

    {
        // reaction 82:  AR2e + ARe => AR + AR2+ + E
        const double k_f = 301107038;
        const double qf = k_f * (sc[39] * sc[42]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[AR_ID] += qdot;
        wdot[ARe_ID] -= qdot;
        wdot[AR2e_ID] -= qdot;
        wdot[AR2p_ID] += qdot;
    }

    {
        // reaction 83:  2 AR2e => 2 AR + AR2+ + E
        const double k_f = 301107038;
        const double qf = k_f * ((sc[42] * sc[42]));
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[E_ID] += qdot;
        wdot[AR_ID] += 2.000000 * qdot;
        wdot[AR2e_ID] -= 2.000000 * qdot;
        wdot[AR2p_ID] += qdot;
    }

    {
        // reaction 85:  AR + ARe => 2 AR
        const double k_f = 1258.62741884;
        const double qf = k_f * (sc[38] * sc[39]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[AR_ID] -= qdot;
        wdot[AR_ID] += 2.000000 * qdot;
        wdot[ARe_ID] -= qdot;
    }

    {
        // reaction 88:  ARe + CO2 => AR + CO + O
        const double k_f = 319173460.28;
        const double qf = k_f * (sc[12] * sc[39]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARe_ID] -= qdot;
    }

    {
        // reaction 89:  ARe + CO => AR + C + O
        const double k_f = 8430997.064;
        const double qf = k_f * (sc[11] * sc[39]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[CO_ID] -= qdot;
        wdot[C_ID] += qdot;
        wdot[AR_ID] += qdot;
        wdot[ARe_ID] -= qdot;
    }

    {
        // reaction 90:  ARe + H2 => AR + 2 H
        const double k_f = 39746129.016;
        const double qf = k_f * (sc[3] * sc[39]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[H2_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARe_ID] -= qdot;
    }

    {
        // reaction 91:  ARe + OH => AR + H + O
        const double k_f = 39746129.016;
        const double qf = k_f * (sc[22] * sc[39]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += qdot;
        wdot[OH_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARe_ID] -= qdot;
    }

    {
        // reaction 92:  ARe + H2O => AR + H + OH
        const double k_f = 126464955.96;
        const double qf = k_f * (sc[24] * sc[39]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[OH_ID] += qdot;
        wdot[H2O_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARe_ID] -= qdot;
    }

    {
        // reaction 93:  ARe + O2 => AR + 2 O
        const double k_f = 126464955.96;
        const double qf = k_f * (sc[9] * sc[39]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += 2.000000 * qdot;
        wdot[O2_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARe_ID] -= qdot;
    }

    {
        // reaction 94:  AR2e + O2 => 2 AR + 2 O
        const double k_f = 27701847.496;
        const double qf = k_f * (sc[9] * sc[42]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += 2.000000 * qdot;
        wdot[O2_ID] -= qdot;
        wdot[AR_ID] += 2.000000 * qdot;
        wdot[AR2e_ID] -= qdot;
    }

    {
        // reaction 95:  ARe => AR
        const double k_f = 0.02;
        const double qf = k_f * (sc[39]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[AR_ID] += qdot;
        wdot[ARe_ID] -= qdot;
    }

    {
        // reaction 96:  AR2e => 2 AR
        const double k_f = 317000;
        const double qf = k_f * (sc[42]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[AR_ID] += 2.000000 * qdot;
        wdot[AR2e_ID] -= qdot;
    }

    {
        // reaction 97:  CO2 + CO2v1 => 2 CO2
        const double k_f = 6443.6906132;
        const double qf = k_f * (sc[12] * sc[13]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2_ID] -= qdot;
        wdot[CO2_ID] += 2.000000 * qdot;
        wdot[CO2v1_ID] -= qdot;
    }

    {
        // reaction 98:  CO + CO2v1 => CO + CO2
        const double k_f = 4504.56128848;
        const double qf = k_f * (sc[11] * sc[13]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO2v1_ID] -= qdot;
    }

    {
        // reaction 99:  CO2v1 + O2 => CO2 + O2
        const double k_f = 4504.56128848;
        const double qf = k_f * (sc[9] * sc[13]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] -= qdot;
        wdot[O2_ID] += qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO2v1_ID] -= qdot;
    }

    {
        // reaction 100:  AR + CO2v1 => AR + CO2
        const double k_f = 4504.56128848;
        const double qf = k_f * (sc[13] * sc[38]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2_ID] += qdot;
        wdot[CO2v1_ID] -= qdot;
        wdot[AR_ID] -= qdot;
        wdot[AR_ID] += qdot;
    }

    {
        // reaction 101:  CO2v1 + H2 => CO2 + H2
        const double k_f = 4504.56128848;
        const double qf = k_f * (sc[3] * sc[13]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO2v1_ID] -= qdot;
    }

    {
        // reaction 102:  CO2 + CO2v2 => 2 CO2
        const double k_f = 5.419926684;
        const double qf = k_f * (sc[12] * sc[14]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2_ID] -= qdot;
        wdot[CO2_ID] += 2.000000 * qdot;
        wdot[CO2v2_ID] -= qdot;
    }

    {
        // reaction 103:  CO + CO2v2 => CO + CO2
        const double k_f = 16.8017727204;
        const double qf = k_f * (sc[11] * sc[14]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO2v2_ID] -= qdot;
    }

    {
        // reaction 104:  CO2v2 + O2 => CO2 + O2
        const double k_f = 16.8017727204;
        const double qf = k_f * (sc[9] * sc[14]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] -= qdot;
        wdot[O2_ID] += qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO2v2_ID] -= qdot;
    }

    {
        // reaction 105:  AR + CO2v2 => AR + CO2
        const double k_f = 16.8017727204;
        const double qf = k_f * (sc[14] * sc[38]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2_ID] += qdot;
        wdot[CO2v2_ID] -= qdot;
        wdot[AR_ID] -= qdot;
        wdot[AR_ID] += qdot;
    }

    {
        // reaction 106:  CO2v2 + H2 => CO2 + H2
        const double k_f = 16.8017727204;
        const double qf = k_f * (sc[3] * sc[14]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO2v2_ID] -= qdot;
    }

    {
        // reaction 107:  CO2 + CO2v2 => CO2 + CO2v1
        const double k_f = 17464.208204;
        const double qf = k_f * (sc[12] * sc[14]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v2_ID] -= qdot;
    }

    {
        // reaction 108:  CO + CO2v2 => CO + CO2v1
        const double k_f = 12224.9457428;
        const double qf = k_f * (sc[11] * sc[14]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v2_ID] -= qdot;
    }

    {
        // reaction 109:  CO2v2 + O2 => CO2v1 + O2
        const double k_f = 12224.9457428;
        const double qf = k_f * (sc[9] * sc[14]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] -= qdot;
        wdot[O2_ID] += qdot;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v2_ID] -= qdot;
    }

    {
        // reaction 110:  AR + CO2v2 => AR + CO2v1
        const double k_f = 12224.9457428;
        const double qf = k_f * (sc[14] * sc[38]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v2_ID] -= qdot;
        wdot[AR_ID] -= qdot;
        wdot[AR_ID] += qdot;
    }

    {
        // reaction 111:  CO2v2 + H2 => CO2v1 + H2
        const double k_f = 12224.9457428;
        const double qf = k_f * (sc[3] * sc[14]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v2_ID] -= qdot;
    }

    {
        // reaction 112:  CO2 + CO2v3 => CO2 + CO2v2
        const double k_f = 464.909266672;
        const double qf = k_f * (sc[12] * sc[15]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO2v2_ID] += qdot;
        wdot[CO2v3_ID] -= qdot;
    }

    {
        // reaction 113:  CO + CO2v3 => CO + CO2v2
        const double k_f = 139.713665632;
        const double qf = k_f * (sc[11] * sc[15]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2v2_ID] += qdot;
        wdot[CO2v3_ID] -= qdot;
    }

    {
        // reaction 114:  CO2v3 + O2 => CO2v2 + O2
        const double k_f = 186.084149484;
        const double qf = k_f * (sc[9] * sc[15]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] -= qdot;
        wdot[O2_ID] += qdot;
        wdot[CO2v2_ID] += qdot;
        wdot[CO2v3_ID] -= qdot;
    }

    {
        // reaction 115:  AR + CO2v3 => AR + CO2v2
        const double k_f = 186.084149484;
        const double qf = k_f * (sc[15] * sc[38]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2v2_ID] += qdot;
        wdot[CO2v3_ID] -= qdot;
        wdot[AR_ID] -= qdot;
        wdot[AR_ID] += qdot;
    }

    {
        // reaction 116:  CO2v3 + H2 => CO2v2 + H2
        const double k_f = 186.084149484;
        const double qf = k_f * (sc[3] * sc[15]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[CO2v2_ID] += qdot;
        wdot[CO2v3_ID] -= qdot;
    }

    {
        // reaction 117:  CO2 + CO2v3 => CO2 + CO2v4
        const double k_f = 3643.3951598;
        const double qf = k_f * (sc[12] * sc[15]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO2v3_ID] -= qdot;
        wdot[CO2v4_ID] += qdot;
    }

    {
        // reaction 118:  CO + CO2v3 => CO + CO2v4
        const double k_f = 1090.00747756;
        const double qf = k_f * (sc[11] * sc[15]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2v3_ID] -= qdot;
        wdot[CO2v4_ID] += qdot;
    }

    {
        // reaction 119:  CO2v3 + O2 => CO2v4 + O2
        const double k_f = 1457.35806392;
        const double qf = k_f * (sc[9] * sc[15]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] -= qdot;
        wdot[O2_ID] += qdot;
        wdot[CO2v3_ID] -= qdot;
        wdot[CO2v4_ID] += qdot;
    }

    {
        // reaction 120:  AR + CO2v3 => AR + CO2v4
        const double k_f = 1457.35806392;
        const double qf = k_f * (sc[15] * sc[38]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2v3_ID] -= qdot;
        wdot[CO2v4_ID] += qdot;
        wdot[AR_ID] -= qdot;
        wdot[AR_ID] += qdot;
    }

    {
        // reaction 121:  CO2v3 + H2 => CO2v4 + H2
        const double k_f = 1457.35806392;
        const double qf = k_f * (sc[3] * sc[15]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[CO2v3_ID] -= qdot;
        wdot[CO2v4_ID] += qdot;
    }

    {
        // reaction 122:  CO2 + CO2v3 => CO2v1 + CO2v2
        const double k_f = 1457.35806392;
        const double qf = k_f * (sc[12] * sc[15]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2_ID] -= qdot;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v2_ID] += qdot;
        wdot[CO2v3_ID] -= qdot;
    }

    {
        // reaction 123:  CO2 + CO2v3 => CO2 + CO2v1
        const double k_f = 1.0237639292;
        const double qf = k_f * (sc[12] * sc[15]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v3_ID] -= qdot;
    }

    {
        // reaction 124:  CO + CO2v3 => CO + CO2v1
        const double k_f = 0.30712917876;
        const double qf = k_f * (sc[11] * sc[15]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v3_ID] -= qdot;
    }

    {
        // reaction 125:  CO2v3 + O2 => CO2v1 + O2
        const double k_f = 0.40950557168;
        const double qf = k_f * (sc[9] * sc[15]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] -= qdot;
        wdot[O2_ID] += qdot;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v3_ID] -= qdot;
    }

    {
        // reaction 126:  AR + CO2v3 => AR + CO2v1
        const double k_f = 0.40950557168;
        const double qf = k_f * (sc[15] * sc[38]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v3_ID] -= qdot;
        wdot[AR_ID] -= qdot;
        wdot[AR_ID] += qdot;
    }

    {
        // reaction 127:  CO2v3 + H2 => CO2v1 + H2
        const double k_f = 0.40950557168;
        const double qf = k_f * (sc[3] * sc[15]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v3_ID] -= qdot;
    }

    {
        // reaction 128:  CO2 + CO2v4 => CO2 + CO2v2
        const double k_f = 26075.8694908;
        const double qf = k_f * (sc[12] * sc[16]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO2v2_ID] += qdot;
        wdot[CO2v4_ID] -= qdot;
    }

    {
        // reaction 129:  CO + CO2v4 => CO + CO2v2
        const double k_f = 18247.0865028;
        const double qf = k_f * (sc[11] * sc[16]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2v2_ID] += qdot;
        wdot[CO2v4_ID] -= qdot;
    }

    {
        // reaction 130:  CO2v4 + O2 => CO2v2 + O2
        const double k_f = 18247.0865028;
        const double qf = k_f * (sc[9] * sc[16]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] -= qdot;
        wdot[O2_ID] += qdot;
        wdot[CO2v2_ID] += qdot;
        wdot[CO2v4_ID] -= qdot;
    }

    {
        // reaction 131:  AR + CO2v4 => AR + CO2v2
        const double k_f = 18247.0865028;
        const double qf = k_f * (sc[16] * sc[38]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2v2_ID] += qdot;
        wdot[CO2v4_ID] -= qdot;
        wdot[AR_ID] -= qdot;
        wdot[AR_ID] += qdot;
    }

    {
        // reaction 132:  CO2v4 + H2 => CO2v2 + H2
        const double k_f = 18247.0865028;
        const double qf = k_f * (sc[3] * sc[16]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[CO2v2_ID] += qdot;
        wdot[CO2v4_ID] -= qdot;
    }

    {
        // reaction 133:  CO2 + CO2v4 => CO2 + CO2v1
        const double k_f = 5.46810381008;
        const double qf = k_f * (sc[12] * sc[16]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v4_ID] -= qdot;
    }

    {
        // reaction 134:  CO + CO2v4 => CO + CO2v1
        const double k_f = 3721.68298968;
        const double qf = k_f * (sc[11] * sc[16]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v4_ID] -= qdot;
    }

    {
        // reaction 135:  CO2v4 + O2 => CO2v1 + O2
        const double k_f = 3721.68298968;
        const double qf = k_f * (sc[9] * sc[16]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] -= qdot;
        wdot[O2_ID] += qdot;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v4_ID] -= qdot;
    }

    {
        // reaction 136:  AR + CO2v4 => AR + CO2v1
        const double k_f = 3721.68298968;
        const double qf = k_f * (sc[16] * sc[38]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v4_ID] -= qdot;
        wdot[AR_ID] -= qdot;
        wdot[AR_ID] += qdot;
    }

    {
        // reaction 137:  CO2v4 + H2 => CO2v1 + H2
        const double k_f = 3721.68298968;
        const double qf = k_f * (sc[3] * sc[16]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[CO2v1_ID] += qdot;
        wdot[CO2v4_ID] -= qdot;
    }

    {
        // reaction 139:  AR+ + O- => AR + O
        const double k_f = 120442815200 * exp((-0.5) * logT);
        const double qf = k_f * (sc[8] * sc[40]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARp_ID] -= qdot;
    }

    {
        // reaction 141:  AR+ + O2- => AR + O2
        const double k_f = 120442815200 * exp((-0.5) * logT);
        const double qf = k_f * (sc[10] * sc[40]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARp_ID] -= qdot;
    }

    {
        // reaction 142:  AR+ + O2- => AR + 2 O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[10] * sc[40]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += 2.000000 * qdot;
        wdot[O2n_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARp_ID] -= qdot;
    }

    {
        // reaction 144:  AR2+ + O- => 2 AR + O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[8] * sc[43]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
        wdot[AR_ID] += 2.000000 * qdot;
        wdot[AR2p_ID] -= qdot;
    }

    {
        // reaction 145:  AR2+ + O2- => 2 AR + O2
        const double k_f = 60221407600;
        const double qf = k_f * (sc[10] * sc[43]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        wdot[AR_ID] += 2.000000 * qdot;
        wdot[AR2p_ID] -= qdot;
    }

    {
        // reaction 146:  AR2+ + O2- => 2 AR + 2 O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[10] * sc[43]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += 2.000000 * qdot;
        wdot[O2n_ID] -= qdot;
        wdot[AR_ID] += 2.000000 * qdot;
        wdot[AR2p_ID] -= qdot;
    }

    {
        // reaction 147:  AR+ + H- => AR + H
        const double k_f = 120442815200 * exp((-0.5) * logT);
        const double qf = k_f * (sc[2] * sc[40]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARp_ID] -= qdot;
    }

    {
        // reaction 149:  AR2+ + H- => 2 AR + H
        const double k_f = 60221407600 * exp((-0.5) * logT);
        const double qf = k_f * (sc[2] * sc[43]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[AR_ID] += 2.000000 * qdot;
        wdot[AR2p_ID] -= qdot;
    }

    {
        // reaction 150:  ARH+ + O- => AR + H + O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[8] * sc[41]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARHp_ID] -= qdot;
    }

    {
        // reaction 152:  ARH+ + O2- => AR + H + O2
        const double k_f = 60221407600;
        const double qf = k_f * (sc[10] * sc[41]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARHp_ID] -= qdot;
    }

    {
        // reaction 153:  ARH+ + O2- => AR + H + 2 O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[10] * sc[41]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += 2.000000 * qdot;
        wdot[O2n_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARHp_ID] -= qdot;
    }

    {
        // reaction 155:  ARH+ + H- => AR + H2
        const double k_f = 3011070380000;
        const double qf = k_f * (sc[2] * sc[41]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[Hn_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[AR_ID] += qdot;
        wdot[ARHp_ID] -= qdot;
    }

    {
        // reaction 156:  ARH+ + H- => AR + 2 H
        const double k_f = 60221407600;
        const double qf = k_f * (sc[2] * sc[41]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[Hn_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARHp_ID] -= qdot;
    }

    {
        // reaction 157:  H2+ + O- => 2 H + O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[5] * sc[8]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[H2p_ID] -= qdot;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
    }

    {
        // reaction 160:  H2+ + O- => H2 + O
        const double k_f = 120442815200 * exp((-0.5) * logT);
        const double qf = k_f * (sc[5] * sc[8]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
    }

    {
        // reaction 161:  H3+ + O- => H + H2 + O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[6] * sc[8]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2_ID] += qdot;
        wdot[H3p_ID] -= qdot;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
    }

    {
        // reaction 162:  O- + OH+ => O + OH
        const double k_f = 120442815200 * exp((-0.5) * logT);
        const double qf = k_f * (sc[8] * sc[23]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[OHp_ID] -= qdot;
    }

    {
        // reaction 163:  O- + OH+ => H + 2 O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[8] * sc[23]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += 2.000000 * qdot;
        wdot[On_ID] -= qdot;
        wdot[OHp_ID] -= qdot;
    }

    {
        // reaction 165:  H2O+ + O- => H2O + O
        const double k_f = 120442815200 * exp((-0.5) * logT);
        const double qf = k_f * (sc[8] * sc[25]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
    }

    {
        // reaction 166:  H2O+ + O- => H + O + OH
        const double k_f = 60221407600;
        const double qf = k_f * (sc[8] * sc[25]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
    }

    {
        // reaction 168:  H3O+ + O- => H + H2O + O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[8] * sc[26]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[H3Op_ID] -= qdot;
    }

    {
        // reaction 169:  H2+ + O2- => H2 + O2
        const double k_f = 120442815200 * exp((-0.5) * logT);
        const double qf = k_f * (sc[5] * sc[10]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
    }

    {
        // reaction 170:  H2+ + O2- => 2 H + O2
        const double k_f = 60221407600;
        const double qf = k_f * (sc[5] * sc[10]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[H2p_ID] -= qdot;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
    }

    {
        // reaction 171:  H2+ + O2- => H2 + 2 O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[5] * sc[10]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[O_ID] += 2.000000 * qdot;
        wdot[O2n_ID] -= qdot;
    }

    {
        // reaction 172:  H2+ + O2- => 2 H + 2 O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[5] * sc[10]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[H2p_ID] -= qdot;
        wdot[O_ID] += 2.000000 * qdot;
        wdot[O2n_ID] -= qdot;
    }

    {
        // reaction 174:  H3+ + O2- => H + H2 + O2
        const double k_f = 60221407600;
        const double qf = k_f * (sc[6] * sc[10]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2_ID] += qdot;
        wdot[H3p_ID] -= qdot;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
    }

    {
        // reaction 175:  O2- + OH+ => O2 + OH
        const double k_f = 120442815200 * exp((-0.5) * logT);
        const double qf = k_f * (sc[10] * sc[23]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[OHp_ID] -= qdot;
    }

    {
        // reaction 176:  O2- + OH+ => H + O + O2
        const double k_f = 60221407600;
        const double qf = k_f * (sc[10] * sc[23]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += qdot;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        wdot[OHp_ID] -= qdot;
    }

    {
        // reaction 177:  O2- + OH+ => 2 O + OH
        const double k_f = 60221407600;
        const double qf = k_f * (sc[10] * sc[23]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += 2.000000 * qdot;
        wdot[O2n_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[OHp_ID] -= qdot;
    }

    {
        // reaction 178:  O2- + OH+ => H + 3 O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[10] * sc[23]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += 3.000000 * qdot;
        wdot[O2n_ID] -= qdot;
        wdot[OHp_ID] -= qdot;
    }

    {
        // reaction 180:  H2O+ + O2- => H2O + O2
        const double k_f = 120442815200 * exp((-0.5) * logT);
        const double qf = k_f * (sc[10] * sc[25]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
    }

    {
        // reaction 181:  H2O+ + O2- => H + O2 + OH
        const double k_f = 60221407600;
        const double qf = k_f * (sc[10] * sc[25]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
    }

    {
        // reaction 182:  H2O+ + O2- => H2O + 2 O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[10] * sc[25]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += 2.000000 * qdot;
        wdot[O2n_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
    }

    {
        // reaction 183:  H2O+ + O2- => H + 2 O + OH
        const double k_f = 60221407600;
        const double qf = k_f * (sc[10] * sc[25]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += 2.000000 * qdot;
        wdot[O2n_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
    }

    {
        // reaction 185:  H3O+ + O2- => H + H2O + O2
        const double k_f = 60221407600;
        const double qf = k_f * (sc[10] * sc[26]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[H3Op_ID] -= qdot;
    }

    {
        // reaction 186:  H3O+ + O2- => H + H2O + 2 O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[10] * sc[26]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += 2.000000 * qdot;
        wdot[O2n_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[H3Op_ID] -= qdot;
    }

    {
        // reaction 187:  CO2+ + O2- => CO + O + O2
        const double k_f = 361328445600;
        const double qf = k_f * (sc[10] * sc[17]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2p_ID] -= qdot;
        //amrex::Print()<<"CO2+ production:"<<wdot[CO2p_ID]<<"\t"<<187<<"\n";
    }

    {
        // reaction 188:  C2O3+ + O2- => CO + CO2 + O2
        const double k_f = 361328445600;
        const double qf = k_f * (sc[10] * sc[18]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2_ID] += qdot;
        wdot[C2O3p_ID] -= qdot;
    }

    {
        // reaction 189:  C2O4+ + O2- => 2 CO2 + O2
        const double k_f = 361328445600;
        const double qf = k_f * (sc[10] * sc[19]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
        wdot[CO2_ID] += 2.000000 * qdot;
        wdot[C2O4p_ID] -= qdot;
    }

    {
        // reaction 190:  H- + H2+ => H + H2
        const double k_f = 120442815200 * exp((-0.5) * logT);
        const double qf = k_f * (sc[2] * sc[5]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[H2p_ID] -= qdot;
    }

    {
        // reaction 191:  H- + H2+ => 3 H
        const double k_f = 60221407600;
        const double qf = k_f * (sc[2] * sc[5]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += 3.000000 * qdot;
        wdot[Hn_ID] -= qdot;
        wdot[H2p_ID] -= qdot;
    }

    {
        // reaction 193:  H- + H3+ => 2 H + H2
        const double k_f = 60221407600;
        const double qf = k_f * (sc[2] * sc[6]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[Hn_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[H3p_ID] -= qdot;
    }

    {
        // reaction 194:  H- + OH+ => 2 H + O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[2] * sc[23]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[Hn_ID] -= qdot;
        wdot[O_ID] += qdot;
        wdot[OHp_ID] -= qdot;
    }

    {
        // reaction 197:  H- + H2O+ => H + H2O
        const double k_f = 120442815200 * exp((-0.5) * logT);
        const double qf = k_f * (sc[2] * sc[25]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
    }

    {
        // reaction 198:  H- + H2O+ => 2 H + OH
        const double k_f = 60221407600;
        const double qf = k_f * (sc[2] * sc[25]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[Hn_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
    }

    {
        // reaction 200:  H- + H3O+ => H2 + H2O
        const double k_f = 138509237480 * exp((-0.5) * logT);
        const double qf = k_f * (sc[2] * sc[26]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[Hn_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[H2O_ID] += qdot;
        wdot[H3Op_ID] -= qdot;
    }

    {
        // reaction 201:  H- + H3O+ => H + H2 + OH
        const double k_f = 138509237480 * exp((-0.5) * logT);
        const double qf = k_f * (sc[2] * sc[26]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[OH_ID] += qdot;
        wdot[H3Op_ID] -= qdot;
    }

    {
        // reaction 204:  CO2+ + H- => CO2 + H
        const double k_f = 116227316668;
        const double qf = k_f * (sc[2] * sc[17]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO2p_ID] -= qdot;
        //amrex::Print()<<"CO2+ production:"<<wdot[CO2p_ID]<<"\t"<<204<<"\n";
    }

    {
        // reaction 205:  CO2+ + H- => CO + H + O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[2] * sc[17]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[Hn_ID] -= qdot;
        wdot[O_ID] += qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2p_ID] -= qdot;
        //amrex::Print()<<"CO2+ production:"<<wdot[CO2p_ID]<<"\t"<<205<<"\n";
    }

    {
        // reaction 206:  CO2+ + CO3- => 2 CO2 + O
        const double k_f = 301107038000;
        const double qf = k_f * (sc[17] * sc[20]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[CO2_ID] += 2.000000 * qdot;
        wdot[CO2p_ID] -= qdot;
        //amrex::Print()<<"CO2+ production:"<<wdot[CO2p_ID]<<"\t"<<206<<"\n";
        wdot[CO3n_ID] -= qdot;
    }

    {
        // reaction 207:  C2O3+ + CO3- => CO + 2 CO2 + O
        const double k_f = 301107038000;
        const double qf = k_f * (sc[18] * sc[20]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2_ID] += 2.000000 * qdot;
        wdot[C2O3p_ID] -= qdot;
        wdot[CO3n_ID] -= qdot;
    }

    {
        // reaction 208:  C2O4+ + CO3- => 3 CO2 + O
        const double k_f = 301107038000;
        const double qf = k_f * (sc[19] * sc[20]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[CO2_ID] += 3.000000 * qdot;
        wdot[C2O4p_ID] -= qdot;
        wdot[CO3n_ID] -= qdot;
    }

    {
        // reaction 209:  CO3- + OH+ => CO2 + O + OH
        const double k_f = 60221407600;
        const double qf = k_f * (sc[20] * sc[23]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO3n_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[OHp_ID] -= qdot;
    }

    {
        // reaction 210:  CO3- + OH+ => CO2 + H + 2 O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[20] * sc[23]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += 2.000000 * qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO3n_ID] -= qdot;
        wdot[OHp_ID] -= qdot;
    }

    {
        // reaction 211:  CO3- + H2O+ => CO2 + H2O + O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[20] * sc[25]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO3n_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
    }

    {
        // reaction 212:  CO3- + H2O+ => CO2 + H + O + OH
        const double k_f = 60221407600;
        const double qf = k_f * (sc[20] * sc[25]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO3n_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
    }

    {
        // reaction 213:  CO3- + H3O+ => CO2 + H + H2O + O
        const double k_f = 60221407600;
        const double qf = k_f * (sc[20] * sc[26]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] += qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO3n_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[H3Op_ID] -= qdot;
    }

    {
        // reaction 214:  AR+ + H2 => ARH+ + H
        const double k_f = 662435483.6;
        const double qf = k_f * (sc[3] * sc[40]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[ARp_ID] -= qdot;
        wdot[ARHp_ID] += qdot;
    }

    {
        // reaction 215:  AR+ + H2 => AR + H2+
        const double k_f = 662435483.6;
        const double qf = k_f * (sc[3] * sc[40]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] -= qdot;
        wdot[H2p_ID] += qdot;
        wdot[AR_ID] += qdot;
        wdot[ARp_ID] -= qdot;
    }

    {
        // reaction 216:  AR+ + H2O => AR + H2O+
        const double k_f = 421549853.2;
        const double qf = k_f * (sc[24] * sc[40]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2O_ID] -= qdot;
        wdot[H2Op_ID] += qdot;
        wdot[AR_ID] += qdot;
        wdot[ARp_ID] -= qdot;
    }

    {
        // reaction 217:  AR+ + H2O => ARH+ + OH
        const double k_f = 180664222.8;
        const double qf = k_f * (sc[24] * sc[40]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[OH_ID] += qdot;
        wdot[H2O_ID] -= qdot;
        wdot[ARp_ID] -= qdot;
        wdot[ARHp_ID] += qdot;
    }

    {
        // reaction 218:  AR+ + CO2 => AR + CO2+
        const double k_f = 457682697.76;
        const double qf = k_f * (sc[12] * sc[40]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2_ID] -= qdot;
        wdot[CO2p_ID] += qdot;
        wdot[AR_ID] += qdot;
        wdot[ARp_ID] -= qdot;
        //amrex::Print()<<"CO2+ production:"<<wdot[CO2p_ID]<<"\t"<<218<<"\n";
    }

    {
        // reaction 219:  AR2+ + H2O => 2 AR + H2O+
        const double k_f = 963542521.6;
        const double qf = k_f * (sc[24] * sc[43]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2O_ID] -= qdot;
        wdot[H2Op_ID] += qdot;
        wdot[AR_ID] += 2.000000 * qdot;
        wdot[AR2p_ID] -= qdot;
    }

    {
        // reaction 220:  AR2+ + H2O => AR + ARH+ + OH
        const double k_f = 240885630.4;
        const double qf = k_f * (sc[24] * sc[43]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[OH_ID] += qdot;
        wdot[H2O_ID] -= qdot;
        wdot[AR_ID] += qdot;
        wdot[ARHp_ID] += qdot;
        wdot[AR2p_ID] -= qdot;
    }

    {
        // reaction 221:  AR2+ + CO2 => 2 AR + CO2+
        const double k_f = 662435483.6;
        const double qf = k_f * (sc[12] * sc[43]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2_ID] -= qdot;
        wdot[CO2p_ID] += qdot;
        wdot[AR_ID] += 2.000000 * qdot;
        wdot[AR2p_ID] -= qdot;
        //amrex::Print()<<"CO2+ production:"<<wdot[CO2p_ID]<<"\t"<<221<<"\n";
    }

    {
        // reaction 222:  ARH+ + H => AR + H2+
        const double k_f = 548014809.16;
        const double qf = k_f * (sc[1] * sc[41]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[H2p_ID] += qdot;
        wdot[AR_ID] += qdot;
        wdot[ARHp_ID] -= qdot;
    }

    {
        // reaction 223:  ARH+ + H2 => AR + H3+
        const double k_f = 301107038;
        const double qf = k_f * (sc[3] * sc[41]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] -= qdot;
        wdot[H3p_ID] += qdot;
        wdot[AR_ID] += qdot;
        wdot[ARHp_ID] -= qdot;
    }

    {
        // reaction 224:  ARH+ + H2O => AR + H3O+
        const double k_f = 1204428152;
        const double qf = k_f * (sc[24] * sc[41]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2O_ID] -= qdot;
        wdot[H3Op_ID] += qdot;
        wdot[AR_ID] += qdot;
        wdot[ARHp_ID] -= qdot;
    }

    {
        // reaction 225:  AR + H2+ => ARH+ + H
        const double k_f = 1385092374.8;
        const double qf = k_f * (sc[5] * sc[38]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[AR_ID] -= qdot;
        wdot[ARHp_ID] += qdot;
    }

    {
        // reaction 226:  AR + H2+ => AR+ + H2
        const double k_f = 132487096.72;
        const double qf = k_f * (sc[5] * sc[38]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[AR_ID] -= qdot;
        wdot[ARp_ID] += qdot;
    }

    {
        // reaction 227:  AR + H3+ => ARH+ + H2
        const double k_f = 6022140.76;
        const double qf = k_f * (sc[6] * sc[38]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H3p_ID] -= qdot;
        wdot[AR_ID] -= qdot;
        wdot[ARHp_ID] += qdot;
    }

    {
        // reaction 228:  H2 + H2+ => H + H3+
        const double k_f = 1264649559.6;
        const double qf = k_f * (sc[3] * sc[5]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2p_ID] -= qdot;
        wdot[H3p_ID] += qdot;
    }

    {
        // reaction 229:  H2+ + O => H + OH+
        const double k_f = 903321114;
        const double qf = k_f * (sc[5] * sc[7]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[O_ID] -= qdot;
        wdot[OHp_ID] += qdot;
    }

    {
        // reaction 230:  H2+ + OH => H2 + OH+
        const double k_f = 457682697.76;
        const double qf = k_f * (sc[5] * sc[22]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[OH_ID] -= qdot;
        wdot[OHp_ID] += qdot;
    }

    {
        // reaction 231:  H2+ + OH => H + H2O+
        const double k_f = 457682697.76;
        const double qf = k_f * (sc[5] * sc[22]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[OH_ID] -= qdot;
        wdot[H2Op_ID] += qdot;
    }

    {
        // reaction 232:  H2+ + H2O => H2 + H2O+
        const double k_f = 2348634896.4;
        const double qf = k_f * (sc[5] * sc[24]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[H2O_ID] -= qdot;
        wdot[H2Op_ID] += qdot;
    }

    {
        // reaction 233:  H2+ + H2O => H + H3O+
        const double k_f = 2047527858.4;
        const double qf = k_f * (sc[5] * sc[24]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[H2O_ID] -= qdot;
        wdot[H3Op_ID] += qdot;
    }

    {
        // reaction 234:  CH4 + H2+ => CH3+ + H + H2
        const double k_f = 1385092374.8;
        const double qf = k_f * (sc[5] * sc[33]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[CH3p_ID] += qdot;
        wdot[CH4_ID] -= qdot;
    }

    {
        // reaction 235:  CH2 + H2+ => CH3+ + H
        const double k_f = 602214076;
        const double qf = k_f * (sc[5] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[CH2_ID] -= qdot;
        wdot[CH3p_ID] += qdot;
    }

    {
        // reaction 236:  CH2 + H2+ => CH2+ + H2
        const double k_f = 602214076;
        const double qf = k_f * (sc[5] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[CH2_ID] -= qdot;
        wdot[CH2p_ID] += qdot;
    }

    {
        // reaction 237:  C + H2+ => CH+ + H
        const double k_f = 1445313782.4;
        const double qf = k_f * (sc[5] * sc[21]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2p_ID] -= qdot;
        wdot[C_ID] -= qdot;
        wdot[CHp_ID] += qdot;
    }

    {
        // reaction 238:  H3+ + O => H + H2O+
        const double k_f = 216797067.36;
        const double qf = k_f * (sc[6] * sc[7]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H3p_ID] -= qdot;
        wdot[O_ID] -= qdot;
        wdot[H2Op_ID] += qdot;
    }

    {
        // reaction 239:  H3+ + O => H2 + OH+
        const double k_f = 505859823.84;
        const double qf = k_f * (sc[6] * sc[7]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H3p_ID] -= qdot;
        wdot[O_ID] -= qdot;
        wdot[OHp_ID] += qdot;
    }

    {
        // reaction 240:  H3+ + OH => H2 + H2O+
        const double k_f = 782878298.8;
        const double qf = k_f * (sc[6] * sc[22]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H3p_ID] -= qdot;
        wdot[OH_ID] -= qdot;
        wdot[H2Op_ID] += qdot;
    }

    {
        // reaction 241:  H2O + H3+ => H2 + H3O+
        const double k_f = 3553063048.4;
        const double qf = k_f * (sc[6] * sc[24]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H3p_ID] -= qdot;
        wdot[H2O_ID] -= qdot;
        wdot[H3Op_ID] += qdot;
    }

    {
        // reaction 242:  CH2 + H3+ => CH3+ + H2
        const double k_f = 1023763929.2;
        const double qf = k_f * (sc[6] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H3p_ID] -= qdot;
        wdot[CH2_ID] -= qdot;
        wdot[CH3p_ID] += qdot;
    }

    {
        // reaction 243:  C + H3+ => CH+ + H2
        const double k_f = 1204428152;
        const double qf = k_f * (sc[6] * sc[21]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[H3p_ID] -= qdot;
        wdot[C_ID] -= qdot;
        wdot[CHp_ID] += qdot;
    }

    {
        // reaction 244:  H2 + OH+ => H + H2O+
        const double k_f = 782878298.8;
        const double qf = k_f * (sc[3] * sc[23]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[OHp_ID] -= qdot;
        wdot[H2Op_ID] += qdot;
    }

    {
        // reaction 245:  OH + OH+ => H2O+ + O
        const double k_f = 421549853.2;
        const double qf = k_f * (sc[22] * sc[23]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[OH_ID] -= qdot;
        wdot[OHp_ID] -= qdot;
        wdot[H2Op_ID] += qdot;
    }

    {
        // reaction 246:  H2O + OH+ => H2O+ + OH
        const double k_f = 963542521.6;
        const double qf = k_f * (sc[23] * sc[24]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[OH_ID] += qdot;
        wdot[OHp_ID] -= qdot;
        wdot[H2O_ID] -= qdot;
        wdot[H2Op_ID] += qdot;
    }

    {
        // reaction 247:  H2O + OH+ => H3O+ + O
        const double k_f = 782878298.8;
        const double qf = k_f * (sc[23] * sc[24]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[OHp_ID] -= qdot;
        wdot[H2O_ID] -= qdot;
        wdot[H3Op_ID] += qdot;
    }

    {
        // reaction 248:  CH4 + OH+ => CH2 + H3O+
        const double k_f = 788900439.56;
        const double qf = k_f * (sc[23] * sc[33]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[OHp_ID] -= qdot;
        wdot[H3Op_ID] += qdot;
        wdot[CH2_ID] += qdot;
        wdot[CH4_ID] -= qdot;
    }

    {
        // reaction 249:  CH2 + OH+ => CH3+ + O
        const double k_f = 289062756.48;
        const double qf = k_f * (sc[23] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[OHp_ID] -= qdot;
        wdot[CH2_ID] -= qdot;
        wdot[CH3p_ID] += qdot;
    }

    {
        // reaction 250:  CH2 + OH+ => CH2+ + OH
        const double k_f = 289062756.48;
        const double qf = k_f * (sc[23] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[OH_ID] += qdot;
        wdot[OHp_ID] -= qdot;
        wdot[CH2_ID] -= qdot;
        wdot[CH2p_ID] += qdot;
    }

    {
        // reaction 251:  C + OH+ => CH+ + O
        const double k_f = 722656891.2;
        const double qf = k_f * (sc[21] * sc[23]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[C_ID] -= qdot;
        wdot[OHp_ID] -= qdot;
        wdot[CHp_ID] += qdot;
    }

    {
        // reaction 252:  H2 + H2O+ => H + H3O+
        const double k_f = 843099706.4;
        const double qf = k_f * (sc[3] * sc[25]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[H2Op_ID] -= qdot;
        wdot[H3Op_ID] += qdot;
    }

    {
        // reaction 253:  H2O+ + OH => H3O+ + O
        const double k_f = 415527712.44;
        const double qf = k_f * (sc[22] * sc[25]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[OH_ID] -= qdot;
        wdot[H2Op_ID] -= qdot;
        wdot[H3Op_ID] += qdot;
    }

    {
        // reaction 254:  CH4 + H2O+ => CH3 + H3O+
        const double k_f = 843099706.4;
        const double qf = k_f * (sc[25] * sc[33]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2Op_ID] -= qdot;
        wdot[H3Op_ID] += qdot;
        wdot[CH3_ID] += qdot;
        wdot[CH4_ID] -= qdot;
    }

    {
        // reaction 255:  CH2 + H2O+ => CH3+ + OH
        const double k_f = 283040615.72;
        const double qf = k_f * (sc[25] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[OH_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
        wdot[CH2_ID] -= qdot;
        wdot[CH3p_ID] += qdot;
    }

    {
        // reaction 256:  CH2 + H2O+ => CH2+ + H2O
        const double k_f = 283040615.72;
        const double qf = k_f * (sc[25] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2O_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
        wdot[CH2_ID] -= qdot;
        wdot[CH2p_ID] += qdot;
    }

    {
        // reaction 257:  C + H2O+ => CH+ + OH
        const double k_f = 662435483.6;
        const double qf = k_f * (sc[21] * sc[25]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[C_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[H2Op_ID] -= qdot;
        wdot[CHp_ID] += qdot;
    }

    {
        // reaction 258:  H2 + H3O+ => H2O + H3+
        const double k_f = 301107038;
        const double qf = k_f * (sc[3] * sc[26]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] -= qdot;
        wdot[H3p_ID] += qdot;
        wdot[H2O_ID] += qdot;
        wdot[H3Op_ID] -= qdot;
    }

    {
        // reaction 259:  CH2 + H3O+ => CH3+ + H2O
        const double k_f = 566081231.44;
        const double qf = k_f * (sc[26] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2O_ID] += qdot;
        wdot[H3Op_ID] -= qdot;
        wdot[CH2_ID] -= qdot;
        wdot[CH3p_ID] += qdot;
    }

    {
        // reaction 260:  CO2+ + H2O => CO2 + H2O+
        const double k_f = 1228516715.04;
        const double qf = k_f * (sc[17] * sc[24]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO2_ID] += qdot;
        wdot[CO2p_ID] -= qdot;
        wdot[H2O_ID] -= qdot;
        wdot[H2Op_ID] += qdot;
        //amrex::Print()<<"CO2+ production:"<<wdot[CO2p_ID]<<"\t"<<260<<"\n";
    }

    {
        // reaction 262:  C2O4+ + CO => C2O3+ + CO2
        const double k_f = 541992668.4;
        const double qf = k_f * (sc[11] * sc[19]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        wdot[C2O3p_ID] += qdot;
        wdot[C2O4p_ID] -= qdot;
    }

    {
        // reaction 265:  CH+ + H2 => CH2+ + H
        const double k_f = 722656891.2;
        const double qf = k_f * (sc[3] * sc[28]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[CHp_ID] -= qdot;
        wdot[CH2p_ID] += qdot;
    }

    {
        // reaction 266:  CH+ + H2O => C + H3O+
        const double k_f = 349284164.08;
        const double qf = k_f * (sc[24] * sc[28]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[C_ID] += qdot;
        wdot[H2O_ID] -= qdot;
        wdot[H3Op_ID] += qdot;
        wdot[CHp_ID] -= qdot;
    }

    {
        // reaction 267:  CH2+ + CH4 => CH3 + CH3+
        const double k_f = 83105542.488;
        const double qf = k_f * (sc[30] * sc[33]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CH2p_ID] -= qdot;
        wdot[CH3_ID] += qdot;
        wdot[CH3p_ID] += qdot;
        wdot[CH4_ID] -= qdot;
    }

    {
        // reaction 268:  CH2+ + H2 => CH3+ + H
        const double k_f = 963542521.6;
        const double qf = k_f * (sc[3] * sc[30]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[CH2p_ID] -= qdot;
        wdot[CH3p_ID] += qdot;
    }

    {
        // reaction 269:  O- + O2 => O + O2-
        const double k_f = 60221407.6;
        const double qf = k_f * (sc[8] * sc[9]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[On_ID] -= qdot;
        wdot[O2_ID] -= qdot;
        wdot[O2n_ID] += qdot;
    }

    {
        // reaction 271:  O + O2- => O- + O2
        const double k_f = 90332111.4 * exp((0.5) * logT);
        const double qf = k_f * (sc[7] * sc[10]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= qdot;
        wdot[On_ID] += qdot;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
    }

    {
        // reaction 272:  H + O2- => H- + O2
        const double k_f = 421549853.2;
        const double qf = k_f * (sc[1] * sc[10]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[Hn_ID] += qdot;
        wdot[O2_ID] += qdot;
        wdot[O2n_ID] -= qdot;
    }

    {
        // reaction 273:  CO3- + O => CO2 + O2-
        const double k_f = 48177126.08;
        const double qf = k_f * (sc[7] * sc[20]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= qdot;
        wdot[O2n_ID] += qdot;
        wdot[CO2_ID] += qdot;
        wdot[CO3n_ID] -= qdot;
    }

    {
        // reaction 282:  H + O2 => O + OH
        const double k_f = 97558680.312 * exp(-(7469.24497053945) * invT);
        const double qf = k_f * (sc[1] * sc[9]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[O_ID] += qdot;
        wdot[O2_ID] -= qdot;
        wdot[OH_ID] += qdot;
    }

    {
        // reaction 290:  H + OH => H2 + O
        const double k_f =
        42154.98532 * exp((2.8) * logT - (1949.96457999329) * invT);
        const double qf = k_f * (sc[1] * sc[22]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[O_ID] += qdot;
        wdot[OH_ID] -= qdot;
    }

    {
        // reaction 296:  H + H2O => H2 + OH
        const double k_f =
        4149254.98364 * exp((1.6) * logT - (9719.12668335236) * invT);
        const double qf = k_f * (sc[1] * sc[24]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[OH_ID] += qdot;
        wdot[H2O_ID] -= qdot;
    }

    {
        // reaction 298:  H2 + O2 => 2 OH
        const double k_f = 1902996480.16 * exp(-(21887.4088791763) * invT);
        const double qf = k_f * (sc[3] * sc[9]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] -= qdot;
        wdot[O2_ID] -= qdot;
        wdot[OH_ID] += 2.000000 * qdot;
    }

    {
        // reaction 299:  H2 + OH => H + H2O
        const double k_f =
        574512.228504 * exp((2) * logT - (1490.02454744777) * invT);
        const double qf = k_f * (sc[3] * sc[22]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[OH_ID] -= qdot;
        wdot[H2O_ID] += qdot;
    }

    {
        // reaction 302:  O + OH => H + O2
        const double k_f =
        10900074.7756 * exp((-0.31) * logT - (-177.132266363261) * invT);
        const double qf = k_f * (sc[7] * sc[22]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] -= qdot;
        wdot[O2_ID] += qdot;
        wdot[OH_ID] -= qdot;
    }

    {
        // reaction 303:  2 OH => H2O + O
        const double k_f =
        33061.5527724 * exp((2.42) * logT - (-969.698515005695) * invT);
        const double qf = k_f * ((sc[22] * sc[22]));
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[OH_ID] -= 2.000000 * qdot;
        wdot[H2O_ID] += qdot;
    }

    {
        // reaction 304:  H2O + O => 2 OH
        const double k_f =
        10056975.0692 * exp((1.14) * logT - (8678.977835134) * invT);
        const double qf = k_f * (sc[7] * sc[24]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= qdot;
        wdot[OH_ID] += 2.000000 * qdot;
        wdot[H2O_ID] -= qdot;
    }

    {
        // reaction 305:  CH2 + CH4 => 2 CH3
        const double k_f = 0.181266436876;
        const double qf = k_f * (sc[29] * sc[33]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CH2_ID] -= qdot;
        wdot[CH3_ID] += 2.000000 * qdot;
        wdot[CH4_ID] -= qdot;
    }

    {
        // reaction 306:  CH4 + H => CH3 + H2
        const double k_f = 177050938.344 * exp(-(6933.31922145743) * invT);
        const double qf = k_f * (sc[1] * sc[33]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[CH3_ID] += qdot;
        wdot[CH4_ID] -= qdot;
    }

    {
        // reaction 307:  CH3 + H => CH2 + H2
        const double k_f = 60221407.6 * exp(-(7599.07487031707) * invT);
        const double qf = k_f * (sc[1] * sc[31]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[CH2_ID] += qdot;
        wdot[CH3_ID] -= qdot;
    }

    {
        // reaction 309:  CH4 + O => CH3 + OH
        const double k_f =
        5010421.11232 * exp((1.56) * logT - (4269.29019268724) * invT);
        const double qf = k_f * (sc[7] * sc[33]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[CH3_ID] += qdot;
        wdot[CH4_ID] -= qdot;
    }

    {
        // reaction 310:  CH3 + O => CH2O + H
        const double k_f = 67447976.512;
        const double qf = k_f * (sc[7] * sc[31]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] -= qdot;
        wdot[CH3_ID] -= qdot;
        wdot[CH2O_ID] += qdot;
    }

    {
        // reaction 311:  CH3 + O => CO + H + H2
        const double k_f = 16861994.128;
        const double qf = k_f * (sc[7] * sc[31]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2_ID] += qdot;
        wdot[O_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[CH3_ID] -= qdot;
    }

    {
        // reaction 312:  CH2 + O => CO + H2
        const double k_f = 33302438.4028;
        const double qf = k_f * (sc[7] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[O_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[CH2_ID] -= qdot;
    }

    {
        // reaction 313:  CH2 + O => CO + 2 H
        const double k_f = 49923546.9004;
        const double qf = k_f * (sc[7] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[O_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[CH2_ID] -= qdot;
    }

    {
        // reaction 314:  CH2 + O2 => CO2 + H2
        const double k_f =
        18006200.8724 * exp((-3.3) * logT - (1439.7028808673) * invT);
        const double qf = k_f * (sc[9] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H2_ID] += qdot;
        wdot[O2_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        wdot[CH2_ID] -= qdot;
    }

    {
        // reaction 315:  CH2 + O2 => CO + H2O
        const double k_f = 855143.98792;
        const double qf = k_f * (sc[9] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[H2O_ID] += qdot;
        wdot[CH2_ID] -= qdot;
    }

    {
        // reaction 316:  CH2 + O2 => CH2O + O
        const double k_f = 324593.386964;
        const double qf = k_f * (sc[9] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[O2_ID] -= qdot;
        wdot[CH2_ID] -= qdot;
        wdot[CH2O_ID] += qdot;
    }

    {
        // reaction 317:  C + O2 => CO + O
        const double k_f = 28304061.572;
        const double qf = k_f * (sc[9] * sc[21]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[O2_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[C_ID] -= qdot;
    }

    {
        // reaction 318:  CH4 + OH <=> CH3 + H2O
        const double k_f =
        81901.114336 * exp((3.04) * logT - (919.880065091028) * invT);
        const double qf = k_f * (sc[22] * sc[33]);
        const double qr =
        k_f * exp(-(g_RT[22] - g_RT[24] - g_RT[31] + g_RT[33])) *
        (sc[24] * sc[31]);
        const double qdot = qf - qr;
        wdot[OH_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[CH3_ID] += qdot;
        wdot[CH4_ID] -= qdot;
    }

    {
        // reaction 679a:  CH4 + HCO => CH2O + CH3
        const double k_f =
        81901.114336 * exp((2.85) * tc[0] - (11300.0) * invT);
        const double qf = k_f * (sc[27] * sc[33]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[HCO_ID] -= qdot;
        wdot[CH3_ID] += qdot;
        wdot[CH4_ID] -= qdot;
        wdot[CH2O_ID] += qdot;
    }

    {
        // reaction 679b:  CH2O + CH3 => CH4 + HCO
        const double k_f = 96.352 * exp((6.1) * tc[0] - (990.0) * invT);
        const double qf = k_f * (sc[31] * sc[34]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[HCO_ID] += qdot;
        wdot[CH3_ID] -= qdot;
        wdot[CH4_ID] += qdot;
        wdot[CH2O_ID] -= qdot;
    }

    {
        // reaction 320:  CH3O + CH4 <=> CH3 + CH3OH
        const double k_f = 157177.873836 * exp(-(4449.44175904533) * invT);
        const double qf = k_f * (sc[33] * sc[35]);
        const double qr =
        k_f * exp(-(-g_RT[31] + g_RT[33] + g_RT[35] - g_RT[36])) *
        (sc[31] * sc[36]);
        const double qdot = qf - qr;
        wdot[CH3_ID] += qdot;
        wdot[CH4_ID] -= qdot;
        wdot[CH3O_ID] -= qdot;
        wdot[CH3OH_ID] += qdot;
    }

    {
        // reaction 321:  CH3 + OH <=> CH2 + H2O
        const double k_f = 72265689.12 * exp(-(1399.94876426873) * invT);
        const double qf = k_f * (sc[22] * sc[31]);
        const double qr =
        k_f * exp(-(g_RT[22] - g_RT[24] - g_RT[29] + g_RT[31])) *
        (sc[24] * sc[29]);
        const double qdot = qf - qr;
        wdot[OH_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[CH2_ID] += qdot;
        wdot[CH3_ID] -= qdot;
    }

    {
        // reaction 322:  CH3 + OH <=> CH2OH + H
        const double k_f =
        927409677.04 * exp((-1.8) * logT - (4059.44884304667) * invT);
        const double qf = k_f * (sc[22] * sc[31]);
        const double qr =
        k_f * exp(-(-g_RT[1] + g_RT[22] + g_RT[31] - g_RT[37])) *
        (sc[1] * sc[37]);
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[OH_ID] -= qdot;
        wdot[CH3_ID] -= qdot;
        wdot[CH2OH_ID] += qdot;
    }

    {
        // reaction 323:  CH3 + OH <=> CH3O + H
        const double k_f =
        1547690.17532 * exp((-0.23) * logT - (7009.30493799394) * invT);
        const double qf = k_f * (sc[22] * sc[31]);
        const double qr =
        k_f * exp(-(-g_RT[1] + g_RT[22] + g_RT[31] - g_RT[35])) *
        (sc[1] * sc[35]);
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[OH_ID] -= qdot;
        wdot[CH3_ID] -= qdot;
        wdot[CH3O_ID] += qdot;
    }

    {
        // reaction 325:  CH3 + HCO => CH4 + CO
        const double k_f = 121045029.276;
        const double qf = k_f * (sc[27] * sc[31]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] += qdot;
        wdot[HCO_ID] -= qdot;
        wdot[CH3_ID] -= qdot;
        wdot[CH4_ID] += qdot;
    }

    {
        // reaction 326:  CH3 + CH3O => CH2O + CH4
        const double k_f = 24088563.04;
        const double qf = k_f * (sc[31] * sc[35]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CH3_ID] -= qdot;
        wdot[CH4_ID] += qdot;
        wdot[CH2O_ID] += qdot;
        wdot[CH3O_ID] -= qdot;
    }

    {
        // reaction 327:  CH2 + CO2 => CH2O + CO
        const double k_f = 23486.348964;
        const double qf = k_f * (sc[12] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[CH2_ID] -= qdot;
        wdot[CH2O_ID] += qdot;
    }

    {
        // reaction 328:  CH2 + OH => CH2O + H
        const double k_f = 18126643.6876;
        const double qf = k_f * (sc[22] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[OH_ID] -= qdot;
        wdot[CH2_ID] -= qdot;
        wdot[CH2O_ID] += qdot;
    }

    {
        // reaction 329:  CH2 + CH2O => CH3 + HCO
        const double k_f = 6022.14076;
        const double qf = k_f * (sc[29] * sc[34]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[HCO_ID] += qdot;
        wdot[CH2_ID] -= qdot;
        wdot[CH3_ID] += qdot;
        wdot[CH2O_ID] -= qdot;
    }

    {
        // reaction 330:  CH2 + HCO => CH3 + CO
        const double k_f = 18126643.6876;
        const double qf = k_f * (sc[27] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] += qdot;
        wdot[HCO_ID] -= qdot;
        wdot[CH2_ID] -= qdot;
        wdot[CH3_ID] += qdot;
    }

    {
        // reaction 331:  CH2 + CH3O => CH2O + CH3
        const double k_f = 18126643.6876;
        const double qf = k_f * (sc[29] * sc[35]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CH2_ID] -= qdot;
        wdot[CH3_ID] += qdot;
        wdot[CH2O_ID] += qdot;
        wdot[CH3O_ID] -= qdot;
    }

    {
        // reaction 700a:  H2 + HCO => CH2O + H
        const double k_f =
        160188.944216 * exp((2) * tc[0] - (8968.83063463752) * invT);
        const double qf = k_f * (sc[3] * sc[27]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[HCO_ID] -= qdot;
        wdot[CH2O_ID] += qdot;
    }

    {
        // reaction 700b:  H + CH2O => H2 + HCO
        const double k_f = 1288708 * exp((1.62) * tc[0] - (1090.0) * invT);
        const double qf = k_f * (sc[1] * sc[34]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[HCO_ID] += qdot;
        wdot[CH2O_ID] -= qdot;
    }

    {
        // reaction 333:  CO2 + H <=> CO + OH
        const double k_f = 151155733.076 * exp(-(13348.325277136) * invT);
        const double qf = k_f * (sc[1] * sc[12]);
        const double qr =
        k_f * exp(-(g_RT[1] - g_RT[11] + g_RT[12] - g_RT[22])) *
        (sc[11] * sc[22]);
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[CO2_ID] -= qdot;
        wdot[OH_ID] += qdot;
    }

    {
        // reaction 335:  H + HCO => CO + H2
        const double k_f = 199935073.232;
        const double qf = k_f * (sc[1] * sc[27]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[CO_ID] += qdot;
        wdot[HCO_ID] -= qdot;
    }

    {
        // reaction 336:  CH3O + H => CH2O + H2
        const double k_f = 13971366.5632;
        const double qf = k_f * (sc[1] * sc[35]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[CH2O_ID] += qdot;
        wdot[CH3O_ID] -= qdot;
    }

    {
        // reaction 337:  CH2O + O => HCO + OH
        const double k_f =
        10719410.5528 * exp((0.57) * logT - (1389.88443095264) * invT);
        const double qf = k_f * (sc[7] * sc[34]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[HCO_ID] += qdot;
        wdot[CH2O_ID] -= qdot;
    }

    {
        // reaction 338:  HCO + O => CO + OH
        const double k_f = 30110703.8;
        const double qf = k_f * (sc[7] * sc[27]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[OH_ID] += qdot;
        wdot[HCO_ID] -= qdot;
    }

    {
        // reaction 339:  HCO + O => CO2 + H
        const double k_f = 30110703.8;
        const double qf = k_f * (sc[7] * sc[27]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        wdot[HCO_ID] -= qdot;
    }

    {
        // reaction 340:  CH3O + O => CH3 + O2
        const double k_f = 21378599.698 * exp(-(239.027916257242) * invT);
        const double qf = k_f * (sc[7] * sc[35]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= qdot;
        wdot[O2_ID] += qdot;
        wdot[CH3_ID] += qdot;
        wdot[CH3O_ID] -= qdot;
    }

    {
        // reaction 341:  CH3O + O => CH2O + OH
        const double k_f = 6022140.76;
        const double qf = k_f * (sc[7] * sc[35]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[CH2O_ID] += qdot;
        wdot[CH3O_ID] -= qdot;
    }

    {
        // reaction 342:  CH3O + CO => CH3 + CO2
        const double k_f = 15717787.3836 * exp(-(5939.46630649311) * invT);
        const double qf = k_f * (sc[11] * sc[35]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        wdot[CH3_ID] += qdot;
        wdot[CH3O_ID] -= qdot;
    }

    {
        // reaction 343:  H2O + HCO => CH2O + OH
        const double k_f =
        514290.820904 * exp((1.35) * logT - (13108.7941442129) * invT);
        const double qf = k_f * (sc[24] * sc[27]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[OH_ID] += qdot;
        wdot[H2O_ID] -= qdot;
        wdot[HCO_ID] -= qdot;
        wdot[CH2O_ID] += qdot;
    }

    {
        // reaction 344:  CH2O + OH => H2O + HCO
        const double k_f =
        2848472.57948 * exp((1.18) * logT - (-224.93784961471) * invT);
        const double qf = k_f * (sc[22] * sc[34]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[OH_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[HCO_ID] += qdot;
        wdot[CH2O_ID] -= qdot;
    }

    {
        // reaction 345:  CH3OH + OH <=> CH3O + H2O
        const double k_f =
        9996753.6616 * exp((1.18) * logT - (853.958681870609) * invT);
        const double qf = k_f * (sc[22] * sc[36]);
        const double qr =
        k_f * exp(-(g_RT[22] - g_RT[24] - g_RT[35] + g_RT[36])) *
        (sc[24] * sc[35]);
        const double qdot = qf - qr;
        wdot[OH_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[CH3O_ID] += qdot;
        wdot[CH3OH_ID] -= qdot;
    }

    {
        // reaction 346:  HCO + OH => CO + H2O
        const double k_f = 101774178.844;
        const double qf = k_f * (sc[22] * sc[27]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] += qdot;
        wdot[OH_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[HCO_ID] -= qdot;
    }

    {
        // reaction 347:  CH3O + OH => CH2O + H2O
        const double k_f = 18126643.6876;
        const double qf = k_f * (sc[22] * sc[35]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[OH_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[CH2O_ID] += qdot;
        wdot[CH3O_ID] -= qdot;
    }

    {
        // reaction 348:  2 HCO => CH2O + CO
        const double k_f = 30110703.8;
        const double qf = k_f * ((sc[27] * sc[27]));
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] += qdot;
        wdot[HCO_ID] -= 2.000000 * qdot;
        wdot[CH2O_ID] += qdot;
    }

    {
        // reaction 349:  CH3O + HCO => CH3OH + CO
        const double k_f = 90332111.4;
        const double qf = k_f * (sc[27] * sc[35]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] += qdot;
        wdot[HCO_ID] -= qdot;
        wdot[CH3O_ID] -= qdot;
        wdot[CH3OH_ID] += qdot;
    }

    {
        // reaction 350:  2 CH3O => CH2O + CH3OH
        const double k_f = 60221407.6;
        const double qf = k_f * ((sc[35] * sc[35]));
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CH2O_ID] += qdot;
        wdot[CH3O_ID] -= 2.000000 * qdot;
        wdot[CH3OH_ID] += qdot;
    }

    {
        // reaction 351:  CH2OH + CH4 <=> CH3 + CH3OH
        const double k_f =
        1011.71964768 * exp((3.1) * logT - (8169.21935267382) * invT);
        const double qf = k_f * (sc[33] * sc[37]);
        const double qr =
        k_f * exp(-(-g_RT[31] + g_RT[33] - g_RT[36] + g_RT[37])) *
        (sc[31] * sc[36]);
        const double qdot = qf - qr;
        wdot[CH3_ID] += qdot;
        wdot[CH4_ID] -= qdot;
        wdot[CH3OH_ID] += qdot;
        wdot[CH2OH_ID] -= qdot;
    }

    {
        // reaction 352:  CH2OH + CH3 => CH2O + CH4
        const double k_f = 2408856.304;
        const double qf = k_f * (sc[31] * sc[37]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CH3_ID] -= qdot;
        wdot[CH4_ID] += qdot;
        wdot[CH2O_ID] += qdot;
        wdot[CH2OH_ID] -= qdot;
    }

    {
        // reaction 353:  CH2 + CH3OH => CH3 + CH3O
        const double k_f =
        674.47976512 * exp((3.1) * logT - (3489.80757735573) * invT);
        const double qf = k_f * (sc[29] * sc[36]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CH2_ID] -= qdot;
        wdot[CH3_ID] += qdot;
        wdot[CH3O_ID] += qdot;
        wdot[CH3OH_ID] -= qdot;
    }

    {
        // reaction 354:  CH2 + CH3OH => CH2OH + CH3
        const double k_f =
        2637.69765288 * exp((3.2) * logT - (3609.57314381725) * invT);
        const double qf = k_f * (sc[29] * sc[36]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CH2_ID] -= qdot;
        wdot[CH3_ID] += qdot;
        wdot[CH3OH_ID] -= qdot;
        wdot[CH2OH_ID] += qdot;
    }

    {
        // reaction 355:  CH2 + CH2OH => CH2O + CH3
        const double k_f = 1210450.29276;
        const double qf = k_f * (sc[29] * sc[37]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CH2_ID] -= qdot;
        wdot[CH3_ID] += qdot;
        wdot[CH2O_ID] += qdot;
        wdot[CH2OH_ID] -= qdot;
    }

    {
        // reaction 356:  CH2OH + H2 <=> CH3OH + H
        const double k_f =
        59980.5219696 * exp((2) * logT - (6726.49717181169) * invT);
        const double qf = k_f * (sc[3] * sc[37]);
        const double qr =
        k_f * exp(-(-g_RT[1] + g_RT[3] - g_RT[36] + g_RT[37])) *
        (sc[1] * sc[36]);
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[H2_ID] -= qdot;
        wdot[CH3OH_ID] += qdot;
        wdot[CH2OH_ID] -= qdot;
    }

    {
        // reaction 357:  CH3OH + H => CH3O + H2
        const double k_f = 39987014.6464 * exp(-(3069.62166140879) * invT);
        const double qf = k_f * (sc[1] * sc[36]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[CH3O_ID] += qdot;
        wdot[CH3OH_ID] -= qdot;
    }

    {
        // reaction 358:  CH2OH + H => CH2O + H2
        const double k_f = 6022140.76;
        const double qf = k_f * (sc[1] * sc[37]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[H2_ID] += qdot;
        wdot[CH2O_ID] += qdot;
        wdot[CH2OH_ID] -= qdot;
    }

    {
        // reaction 359:  CH2OH + H => CH3O + H
        const double k_f = 174642082.04 * exp((0.04) * logT);
        const double qf = k_f * (sc[1] * sc[37]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[H_ID] += qdot;
        wdot[CH3O_ID] += qdot;
        wdot[CH2OH_ID] -= qdot;
    }

    {
        // reaction 360:  CH3O + H => CH3OH
        const double k_f =
        95752038.084 * exp((0.24) * logT - (-3181.83897788324) * invT);
        const double qf = k_f * (sc[1] * sc[35]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] -= qdot;
        wdot[CH3O_ID] -= qdot;
        wdot[CH3OH_ID] += qdot;
    }

    {
        // reaction 361:  CH3OH + O => CH2OH + OH
        const double k_f = 4281742.08036 * exp(-(2026.95672986141) * invT);
        const double qf = k_f * (sc[7] * sc[36]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[CH3OH_ID] -= qdot;
        wdot[CH2OH_ID] += qdot;
    }

    {
        // reaction 362:  CH3OH + O => CH3O + OH
        const double k_f = 9996753.6616 * exp(-(2359.58294595833) * invT);
        const double qf = k_f * (sc[7] * sc[36]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[CH3O_ID] += qdot;
        wdot[CH3OH_ID] -= qdot;
    }

    {
        // reaction 363:  CH2OH + O => CH2O + OH
        const double k_f = 42154985.32;
        const double qf = k_f * (sc[7] * sc[37]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[CH2O_ID] += qdot;
        wdot[CH2OH_ID] -= qdot;
    }

    {
        // reaction 364:  CH3OH + OH => CH2OH + H2O
        const double k_f =
        207161.642144 * exp((2.8) * logT - (-209.841349640568) * invT);
        const double qf = k_f * (sc[22] * sc[36]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[OH_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[CH3OH_ID] -= qdot;
        wdot[CH2OH_ID] += qdot;
    }

    {
        // reaction 365:  CH2OH + OH => CH2O + H2O
        const double k_f = 24088563.04;
        const double qf = k_f * (sc[22] * sc[37]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[OH_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[CH2O_ID] += qdot;
        wdot[CH2OH_ID] -= qdot;
    }

    {
        // reaction 366:  CH2O + CH2OH <=> CH3OH + HCO
        const double k_f =
        46490.9266672 * exp((2.8) * logT - (2949.85609494727) * invT);
        const double qf = k_f * (sc[34] * sc[37]);
        const double qr =
        k_f * exp(-(-g_RT[27] + g_RT[34] - g_RT[36] + g_RT[37])) *
        (sc[27] * sc[36]);
        const double qdot = qf - qr;
        wdot[HCO_ID] += qdot;
        wdot[CH2O_ID] -= qdot;
        wdot[CH3OH_ID] += qdot;
        wdot[CH2OH_ID] -= qdot;
    }

    {
        // reaction 367:  CH2OH + HCO => 2 CH2O
        const double k_f = 181266436.876;
        const double qf = k_f * (sc[27] * sc[37]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[HCO_ID] -= qdot;
        wdot[CH2O_ID] += 2.000000 * qdot;
        wdot[CH2OH_ID] -= qdot;
    }

    {
        // reaction 368:  CH2OH + HCO => CH3OH + CO
        const double k_f = 121045029.276;
        const double qf = k_f * (sc[27] * sc[37]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] += qdot;
        wdot[HCO_ID] -= qdot;
        wdot[CH3OH_ID] += qdot;
        wdot[CH2OH_ID] -= qdot;
    }

    {
        // reaction 369:  CH3O + CH3OH <=> CH2OH + CH3OH
        const double k_f = 301107.038 * exp(-(2049.60147982262) * invT);
        const double qf = k_f * (sc[35] * sc[36]);
        const double qr =
        k_f * exp(-(g_RT[35] + g_RT[36] - g_RT[36] - g_RT[37])) *
        (sc[36] * sc[37]);
        const double qdot = qf - qr;
        wdot[CH3O_ID] -= qdot;
        wdot[CH3OH_ID] -= qdot;
        wdot[CH3OH_ID] += qdot;
        wdot[CH2OH_ID] += qdot;
    }

    {
        // reaction 370:  CH2OH + CH3O => CH2O + CH3OH
        const double k_f = 24088563.04;
        const double qf = k_f * (sc[35] * sc[37]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CH2O_ID] += qdot;
        wdot[CH3O_ID] -= qdot;
        wdot[CH3OH_ID] += qdot;
        wdot[CH2OH_ID] -= qdot;
    }

    {
        // reaction 371:  2 CH2OH => CH2O + CH3OH
        const double k_f = 4817712.608;
        const double qf = k_f * ((sc[37] * sc[37]));
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CH2O_ID] += qdot;
        wdot[CH3OH_ID] += qdot;
        wdot[CH2OH_ID] -= 2.000000 * qdot;
    }

    {
        // reaction 372:  CH2 + O => H + HCO
        const double k_f = 30170925.2076;
        const double qf = k_f * (sc[7] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] -= qdot;
        wdot[HCO_ID] += qdot;
        wdot[CH2_ID] -= qdot;
    }

    {
        // reaction 373:  CH3 + O => CH3O
        const double k_f =
        73470117.272 * exp((0.05) * logT - (-68.4374665494419) * invT);
        const double qf = k_f * (sc[7] * sc[31]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] -= qdot;
        wdot[CH3_ID] -= qdot;
        wdot[CH3O_ID] += qdot;
    }

    {
        // reaction 375:  C + CO2 => 2 CO
        const double k_f = 602.214076;
        const double qf = k_f * (sc[12] * sc[21]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[CO_ID] += 2.000000 * qdot;
        wdot[CO2_ID] -= qdot;
        wdot[C_ID] -= qdot;
    }

    {
        // reaction 377:  CO + O2 => CO2 + O
        const double k_f = 2529299.1192 * exp(-(2399.84027922271) * invT);
        const double qf = k_f * (sc[9] * sc[11]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O_ID] += qdot;
        wdot[O2_ID] -= qdot;
        wdot[CO_ID] -= qdot;
        wdot[CO2_ID] += qdot;
    }

    {
        // reaction 379:  CH2 + O2 => CO2 + 2 H
        const double k_f = 6503912.0208 * exp(-(757.844298701908) * invT);
        const double qf = k_f * (sc[9] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += 2.000000 * qdot;
        wdot[O2_ID] -= qdot;
        wdot[CO2_ID] += qdot;
        wdot[CH2_ID] -= qdot;
    }

    {
        // reaction 380:  CH2 + O2 => CO + H + OH
        const double k_f = 6503912.0208 * exp(-(757.844298701908) * invT);
        const double qf = k_f * (sc[9] * sc[29]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O2_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[OH_ID] += qdot;
        wdot[CH2_ID] -= qdot;
    }

    {
        // reaction 381:  CH3OH + OH => CH2O + H + H2O
        const double k_f =
        662435.4836 * exp((1.44) * logT - (6799.46358835337) * invT);
        const double qf = k_f * (sc[22] * sc[36]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[OH_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[CH2O_ID] += qdot;
        wdot[CH3OH_ID] -= qdot;
    }

    {
        // reaction 382:  CH3 + O2 => CH2O + OH
        const double k_f = 340250.95294 * exp(-(4499.7634256258) * invT);
        const double qf = k_f * (sc[9] * sc[31]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] -= qdot;
        wdot[OH_ID] += qdot;
        wdot[CH3_ID] -= qdot;
        wdot[CH2O_ID] += qdot;
    }

    {
        // reaction 383:  CH3 + O2 => H2O + HCO
        const double k_f = 999675.36616;
        const double qf = k_f * (sc[9] * sc[31]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[O2_ID] -= qdot;
        wdot[H2O_ID] += qdot;
        wdot[HCO_ID] += qdot;
        wdot[CH3_ID] -= qdot;
    }

    {
        // reaction 384:  CH2O + O => CO + H + OH
        const double k_f = 60221407.6;
        const double qf = k_f * (sc[7] * sc[34]);
        const double qr = 0.0;
        const double qdot = qf - qr;
        wdot[H_ID] += qdot;
        wdot[O_ID] -= qdot;
        wdot[CO_ID] += qdot;
        wdot[OH_ID] += qdot;
        wdot[CH2O_ID] -= qdot;
    }

    double totalmass=0.0;
    double totalchrg=0.0;
    int chrgspec[NUM_SPECIES] = {0};
    CKCHRG(chrgspec);
    for (int i = 0; i < 44; ++i)
    {
        totalmass += h_global_mw[i]*wdot[i];
        totalchrg += chrgspec[i]*wdot[i];
    }
    double tol=1e-2;
    int problem=0;
    if(amrex::Math::abs(totalmass)>tol)
    {
        amrex::Print()<<"total mass error:"<<totalmass<<"\n";
        problem=1;
    }
    if(amrex::Math::abs(totalchrg)>tol)
    {
        amrex::Print()<<"total charge error:"<<totalchrg<<"\n";
        problem=1;
    }
    //carbon balance
    int elementcomp[NUM_GAS_SPECIES*5]={0};
    amrex::Vector<std::string> specnames;
    CKNCF(elementcomp);
    CKSYMS_STR(specnames);
    double totalcprod=0.0;
    double totalhprod=0.0;
    double totaloprod=0.0;
    for(int i=0;i<NUM_GAS_SPECIES;i++)
    {
           totalcprod += elementcomp[i*5+2]*wdot[i];
           totalhprod += elementcomp[i*5+4]*wdot[i];
           totaloprod += elementcomp[i*5+1]*wdot[i];
    }
    if(amrex::Math::abs(totalcprod)>tol)
    {
        amrex::Print()<<"total carbon error:"<<totalcprod<<"\n";
        problem=1;
    }
    if(problem)
    {
        amrex::Print()<<"mass/charge/C/H/O error:"<<totalmass<<"\t"<<
        totalchrg<<"\t"<<totalcprod<<"\t"<<totalhprod<<"\t"<<totaloprod<<"\n";
        amrex::Abort("mass/charge/carbon error");
    }

}
