{
    if (reinit_int > 0 && (step % reinit_int == 0)) {
        amrex::Print() << "Reinitializing at step " << step << "\n";
        reinitialize_phi(phi, geom, reinit_iters, reinit_dtau);
    }

    if (plot_int > 0 && (step % plot_int == 0)) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "plt%04d", step);
        Vector<std::string> names = {"phi"};
        WriteSingleLevelPlotfile(buf, phi, names, geom, step*dt, step);
        amrex::Print() << "Wrote " << buf << "\n";
    }   
}