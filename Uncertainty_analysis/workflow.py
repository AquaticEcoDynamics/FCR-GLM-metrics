import os
import sys
import shutil
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0,".")
import pyemu


def reweight_and_refosm(org_d="model_files"):
    """adjust the weights and rerun the fosm calcs
    """

    # load the control file
    pst = pyemu.Pst(os.path.join(org_d,"glm3.pst"))
    # calculate the squared weighted residuals using the weights and obsvals in the control
    # file and the modelled values in the res file
    pst.res.loc[:,"org_swr"] = ((pst.observation_data.obsval.values - pst.res.loc[pst.obs_names,"modelled"]) * pst.observation_data.weight.values)**2

    # check that these numbers agree
    assert np.abs(pst.phi - pst.res.org_swr.sum()) < 1.0e-5

    # copy the original obs data for later
    org_obs = pst.observation_data.copy()

    # load the prior par cov matrix from the unc file
    parcov = pyemu.Cov.from_uncfile(os.path.join(org_d,"glm3_prior.unc"))

    #create a schur object
    org_schur = pyemu.Schur(jco=os.path.join(org_d,"glm3.jco"),pst=pst,parcov=parcov)

    #run the predunc calcs and get a parameter summary dataframe
    org_parsum = org_schur.get_parameter_summary()

    # reset the existing weights to be measurement based
    obs =  pst.observation_data
    obs.loc[:,"standard_deviation"] = np.nan

    # the stdev of the temp measurements
    temp_data_std = 0.52

    obs.loc[obs.apply(lambda x: x.obgnme.startswith("temp_") and x.weight > 0,axis=1),"weight"] = 1.0/temp_data_std
    obs.loc[obs.apply(lambda x: x.obgnme.startswith("temp_") and x.weight > 0,axis=1),"standard_deviation"] = temp_data_std

    # the stdev of the oxy measurements
    oxy_data_std = 16.43

    obs.loc[obs.apply(lambda x: x.obgnme.startswith("oxy_") and x.weight > 0,axis=1),"weight"] = 1.0/oxy_data_std
    obs.loc[obs.apply(lambda x: x.obgnme.startswith("oxy_") and x.weight > 0,axis=1),"standard_deviation"] = oxy_data_std

    # the stdev of the mom measurements
    mom_data_std = 17.6

    obs.loc[obs.apply(lambda x: x.obgnme.startswith("mom") and x.weight > 0,axis=1),"weight"] = 1.0/mom_data_std
    obs.loc[obs.apply(lambda x: x.obgnme.startswith("mom") and x.weight > 0,axis=1),"standard_deviation"] = mom_data_std

    # the stdev of the af measurements
    af_data_std = 4.37

    obs.loc[obs.apply(lambda x: x.obgnme.startswith("af") and x.weight > 0,axis=1),"weight"] = 1.0/af_data_std
    obs.loc[obs.apply(lambda x: x.obgnme.startswith("af") and x.weight > 0,axis=1),"standard_deviation"] = af_data_std

    # the stdev of the td measurements
    td_data_std = 0.35

    obs.loc[obs.apply(lambda x: x.obgnme.startswith("thermo") and x.weight > 0,axis=1),"weight"] = 1.0/td_data_std
    obs.loc[obs.apply(lambda x: x.obgnme.startswith("thermo") and x.weight > 0,axis=1),"standard_deviation"] = td_data_std

    # the stdev of the ss measurements
    ss_data_std = 2.89

    obs.loc[obs.apply(lambda x: x.obgnme.startswith("schmidt") and x.weight > 0,axis=1),"weight"] = 1.0/ss_data_std
    obs.loc[obs.apply(lambda x: x.obgnme.startswith("schmidt") and x.weight > 0,axis=1),"standard_deviation"] = ss_data_std


    # the metrics...
    #obs.loc[obs.apply(lambda x: not x.obgnme.startswith("temp_") and not x.obgnme.startswith("oxy_"),axis=1),"weight"] = 1.0

    pst.write(os.path.join(org_d,"glm3_reweight.pst"),version=2)

    org_phi = pst.phi

    # now do the discrepancy weight adjustment
    pst.adjust_weights_discrepancy()
    print(org_phi,pst.phi)

    # now redo the schur calculations
    new_schur = pyemu.Schur(jco=os.path.join(org_d,"glm3.jco"),pst=pst,parcov=parcov)
    new_parsum = new_schur.get_parameter_summary()

    fig,axes = plt.subplots(3,1,figsize=(20,8))
    org_parsum.loc[:,["prior_var","post_var"]].plot(kind='bar',ax=axes[0])
    axes[0].set_title("original prior vs posterior variance",loc="left")
    new_parsum.loc[:,["prior_var","post_var"]].plot(kind='bar',ax=axes[1])
    axes[1].set_title("new prior vs posterior variance",loc="left")
    s = new_schur.xtqx.s.x
    eigthresh_value = s / s.max()
    xvals = np.arange(pst.npar_adj)+1
    axes[2].plot(xvals,np.log10(eigthresh_value),label="singular value ratio to max")
    axes[2].set_xlabel("singular component")
    axes[2].set_ylabel("$log_{10}$ singular value/max(singular value)")
    xlim = axes[2].get_xlim()
    axes[2].plot(xlim,[-7,-7],"k--",lw=2,label="typical truncation level")
    axes[2].set_xlim(xlim)
    axes[2].set_title("singular spectrum ratio of xtqx",loc="left")
    axes[2].set_xticks(xvals)
    axes[2].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(org_d,"sum.pdf"))


    # now run the identifiability calcs knowing that the number
    # of solution space components is npar - 1
    ev = pyemu.ErrVar(jco=os.path.join(org_d,"glm3.jco"),pst=pst,parcov=parcov)
    id_df = ev.get_identifiability_dataframe(singular_value=pst.npar_adj-1)
    fig,ax = plt.subplots(1,1,figsize=(10,3))
    id_df.ident.plot(kind="bar",ax=ax)
    ax.set_title("identifiability",loc="left")
    plt.tight_layout()
    plt.savefig(os.path.join(org_d,"ident.pdf"))
    plt.close(fig)


    # run the schur calcs using each obs group as if it is the only obs...
    from matplotlib.backends.backend_pdf import PdfPages
    obs_bak = pst.observation_data.copy()
    with PdfPages(os.path.join(org_d,"obsgroup.pdf")) as pdf:
        for g in pst.nnz_obs_groups:
            obs = obs_bak.copy()
            keep = obs.loc[obs.apply(lambda x: x.obgnme==g and x.weight > 0,axis=1),"obsnme"]
            obs.loc[:,"weight"] = 0
            obs.loc[keep,"weight"] = org_obs.loc[keep,"weight"]
            new_schur = pyemu.Schur(jco=os.path.join(org_d,"glm3.jco"),pst=pst,parcov=parcov)
            new_parsum = new_schur.get_parameter_summary()
            fig,axes = plt.subplots(2,1,figsize=(20,5))
            new_parsum.loc[:,["prior_var","post_var"]].plot(kind='bar',ax=axes[0])
            axes[0].set_title("new prior vs posterior variance with only obs group {0}".format(g),loc="left")
            s = new_schur.xtqx.s.x
            eigthresh_value = s / s.max()
            xvals = np.arange(pst.npar_adj)+1
            axes[1].plot(xvals,np.log10(eigthresh_value),label="singular value ratio to max")
            axes[1].set_xlabel("singular component")
            axes[1].set_ylabel("$log_{10}$ singular value/max(singular value)")
            xlim = axes[1].get_xlim()
            axes[1].plot(xlim,[-7,-7],"k--",lw=2,label="typical truncation level")
            axes[1].set_xlim(xlim)
            axes[1].set_title("singular spectrum ratio of xtqx",loc="left")
            axes[1].set_xticks(xvals)
            axes[1].legend()
            plt.tight_layout()
            pdf.savefig()
            plt.close(fig)
            print(g)


def prep_for_ies(org_d="model_files"):
    """prep for an ies run
    """

    # the location of the ies prep'd model+pest files
    t_d = "template"
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(org_d,t_d)

    # the number of realizations to use
    num_reals = 300

    # load the reweighted pst
    pst_name = os.path.join(t_d,"glm3_reweight.pst")
    if not os.path.exists(pst_name):
        raise Exception("need to run 'reweight_and_refosm()' first")
    pst = pyemu.Pst(pst_name)

    # set some extra data on the obs dataframe
    obs = pst.observation_data
    obs.loc[:,"day"] = obs.obsnme.apply(lambda x: int(x.split("_")[2]))

    #these upper and lower bounds are important for when we
    #draw noise realizations to make sure we have
    #physically plausible values in the noise realizations
    obs.loc[:,"lower_bound"] = np.nan
    obs.loc[:,"upper_bound"] = np.nan
    obs.loc[obs.obsnme.str.contains("ox"),"lower_bound"] = 0.0
    obs.loc[obs.obsnme.str.contains("af"),"lower_bound"] = 0.0
    obs.loc[obs.obsnme.str.contains("ss"),"lower_bound"] = 0.0
    obs.loc[obs.obsnme.str.contains("_td_"),"lower_bound"] = 0.0

    # prep for autocorrelated noise draws
    groups = obs.obgnme.unique()
    groups.sort()
    sd = {}

    for g in groups:
        v = pyemu.geostats.ExpVario(contribution=1.0,a=60) # a 60 day correlation range?
        gs = pyemu.geostats.GeoStruct(variograms=v,name=g)
        gobs = obs.loc[obs.obgnme==g,:].copy()
        gobs.sort_values(by="day",inplace=True)
        sd[gs] = gobs.obsnme.tolist()
    # draw the autocorrelated noise realizations
    oe = pyemu.helpers.autocorrelated_draw(pst,struct_dict=sd,time_distance_col="day",num_reals=300)

    #enforce bounds in the realized values
    for oname,lb,ub in zip(obs.obsnme,obs.lower_bound,obs.upper_bound):
        if ~np.isnan(lb):
            oe.loc[:,oname] = oe.loc[:,oname].apply(lambda x: max(lb,x))
        if ~np.isnan(ub):
            oe.loc[:,oname] = oe.loc[:,oname].apply(lambda x: min(ub,x))
    # save
    oe.to_csv(os.path.join(t_d,"noise.csv"))
    pst.pestpp_options["ies_obs_en"] = "noise.csv"

    #write phi factors file that tells ies to rebalance the
    # weights on-the-fly for us.
    # I assumed 30% temp, 30% oxy, and 10% for each metrix type...
    with open(os.path.join(t_d,"phi_facs.dat"),'w') as f:
        f.write("temp 0.3\n")
        f.write("oxy 0.3\n")
        f.write("thermo 0.1\n")
        f.write("schmidt 0.1\n")
        f.write("mom 0.1\n")
        f.write("af 0.1\n")
    #some ies options...
    pst.pestpp_options["ies_phi_factor_file"] = "phi_facs.dat"
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_subset_size"] = -10
    pst.pestpp_options["ies_bad_phi_sigma"] = 2.0
    pst.control_data.noptmax = 3

    #check that ies likes the inputs
    pst.pestpp_options["debug_parse_only"] = True
    pst_name = pst_name.replace(".pst","_ies.pst")
    pst.write(pst_name,version=2)
    pyemu.os_utils.run("pestpp-ies {0}".format(os.path.split(pst_name)[1]),cwd=t_d)

    # now do the ever important test run
    pst.pestpp_options["debug_parse_only"] = False
    pst.control_data.noptmax = -2
    pst.write(pst_name,version=2)
    pyemu.os_utils.run("pestpp-ies {0}".format(os.path.split(pst_name)[1]),cwd=t_d)

    #now save for actually ies run
    pst.control_data.noptmax = 3
    pst.write(pst_name,version=2)

    # plot the noise realizations
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(os.path.join(t_d,"noise.pdf")) as pdf:
        for g in groups:
            gobs = obs.loc[obs.obgnme==g,:].copy()
            gobs.sort_values(by="day",inplace=True)
            fig,ax = plt.subplots(1,1,figsize=(8,3))
            days = gobs.day.values
            names = gobs.obsnme.values
            [ax.plot(days,oe.loc[idx,names],"r",lw=0.1,alpha=0.5) for idx in oe.index]
            ax.set_title(g,loc="left")
            plt.tight_layout()
            pdf.savefig()
            plt.close(fig)


def run_ies_locally(t_d="template",m_d="master_ies",num_workers=10):
    """run pestpp-ies locally in parallel

    Note: set the number of workers to appropriate for what your machine can handle!
    """
    pyemu.os_utils.start_workers(t_d,"pestpp-ies","glm3_reweight_ies.pst",num_workers=num_workers,
                                 worker_root=".",master_dir=m_d)

if __name__ == "__main__":
    reweight_and_refosm()
    prep_for_ies()
    run_ies_locally()
