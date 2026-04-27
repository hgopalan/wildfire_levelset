{
    if (inputs.reinit_int > 0 && (step % inputs.reinit_int == 0)) {
        amrex::Print() << "Reinitializing at step " << step << "\n";
        reinitialize_phi(phi, geom, inputs.reinit_iters, inputs.reinit_dtau);
    }

    if (inputs.plot_int > 0 && (step % inputs.plot_int == 0)) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "plt%04d", step);
        Vector<std::string> names = {"phi"};
        WriteSingleLevelPlotfile(buf, phi, names, geom, step*dt, step);
        amrex::Print() << "Wrote " << buf << "\n";
    }   
}