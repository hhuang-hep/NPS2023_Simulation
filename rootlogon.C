{
    char *inc_nps_soft = gSystem->ExpandPathName("$folder_NPS_SOFT");
    gInterpreter->AddIncludePath(inc_nps_soft);
    delete [] inc_nps_soft;
}
