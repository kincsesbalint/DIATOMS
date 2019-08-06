def stdventr2diff_workflow(SinkTag="std2diff", wf_name="HOXventr2diff"):


    """
        Part of DIATOMS project. The registration of T1 images and MNI to evaluate a coarse lateral ventricle mask.

        Workflow inputs:
            :param highres_brain: The oriented highres brain extracted T1w image.
            :param highres: The oriented high res T1w image.
            :param reference_brain: the standard used as reference (default MNI152_1mm)
            :param reference_ventricle: the standard ventricle mask, made from the Harvard Oxford Lateral Ventricle mask (1mm).
            :param highres2diff: Highres to diffusion affine matrix.
            :param bbr_schedule: Parameters which specifies BBR options.
            :param SinkDir:
            :param SinkTag: The output directory in which the returned images (see workflow outputs) could be found.

        Workflow outputs:




            :return: stdventr2diff_workflow - workflow
                func="/home/balint/Dokumentumok/phd/essen/PAINTER/probe/s002/func_data.nii.gz",
                 skull="/home/balint/Dokumentumok/phd/essen/PAINTER/probe/MS001/highres.nii.gz",
                 anat_wm_segmentation="/home/balint/Dokumentumok/phd/essen/PAINTER/probe/anat_preproc/fast/fast__prob_2.nii.gz",



        Balint Kincses
        kincses.balint@med.u-szeged.hu
        2019


        """
    import os
    import nipype.pipeline as pe
    from nipype.interfaces.utility import Function
    import nipype.interfaces.utility as utility
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.io as io

    SinkDir = os.path.abspath(globals._SinkDir_ + "/" + SinkTag)
    if not os.path.exists(SinkDir):
        os.makedirs(SinkDir)

    # Define inputs of the workflow
    inputspec = pe.Node(utility.IdentityInterface(fields=['highres_brain',
                                                       'highres',
                                                       'reference_brain',
                                                       'reference_ventricle',
                                                       'fnirt_config',
                                                       'highres2diff']),
                        name='inputspec')

    inputspec.inputs.reference_brain = globals._FSLDIR_ + globals._brainref
    inputspec.inputs.fnirt_config = "T1_2_MNI152_2mm"

    # Linear registration node
    linear_reg = pe.MapNode(interface=fsl.FLIRT(),
                         iterfield=['in_file'],
                         name='linear_reg_0')
    linear_reg.inputs.cost = 'corratio'

    # Non-linear registration node
    nonlinear_reg = pe.MapNode(interface=fsl.FNIRT(),
                               iterfield=['in_file', 'affine_file'],
                               name='nonlinear_reg_1')
    nonlinear_reg.inputs.fieldcoeff_file = True
    nonlinear_reg.inputs.jacobian_file = True

    # Calculate inverse of the nonlinear warping field
    inv_fnirt_xfm = pe.MapNode(interface=fsl.utils.InvWarp(),
                               iterfield=['warp', 'reference'],
                               name="inv_nonlinear_xfm")

    # Applying inverted warp field
    brain_warp = pe.MapNode(interface=fsl.ApplyWarp(),
                            iterfield=['in_file', 'field_file','postmat'],
                            name='brain_warp')

    # Threshold and binarise the converted ventricle mask in diffusion space
    thresholdbin = pe.MapNode(interface=fsl.Threshold(),
                            iterfield=['in_file', 'thresh','args'],
                            name='thresholdbin')
    thresholdbin.inputs.args='-bin'
    thresholdbin.inputs.thresh=0.5

    # Save outputs which are important
    ds = pe.Node(interface=io.DataSink(), name='ds')
    ds.inputs.base_directory = SinkDir
    ds.inputs.regexp_substitutions = [("(\/)[^\/]*$", ".nii.gz")]

    # Define outputs of the workflow
    outputspec = pe.Node(utility.IdentityInterface(fields=['output_ventricle']),
                         name='outputspec')

    # Create workflow nad connect nodes
    analysisflow = pe.Workflow(name=wf_name)
    analysisflow.connect(inputspec, 'highres_brain', linear_reg, 'in_file')
    analysisflow.connect(inputspec, 'reference_brain', linear_reg, 'reference')

    analysisflow.connect(inputspec, 'highres_brain', nonlinear_reg, 'in_file')
    analysisflow.connect(inputspec, 'reference_brain', nonlinear_reg, 'ref_file')
    # FNIRT parameters are specified by FSL config file
    # ${FSLDIR}/etc/flirtsch/TI_2_MNI152_2mm.cnf (or user-specified)
    analysisflow.connect(inputspec, 'fnirt_config', nonlinear_reg, 'config_file')
    analysisflow.connect(linear_reg, 'out_matrix_file', nonlinear_reg, 'affine_file')

    analysisflow.connect(nonlinear_reg, 'fieldcoeff_file', inv_fnirt_xfm, 'warp')
    analysisflow.connect(inputspec, 'reference_brain', inv_fnirt_xfm, 'reference')

    analysisflow.connect(inputspec, 'reference_ventricle', brain_warp, 'in_file')
    analysisflow.connect(inputspec, 'reference_brain', brain_warp, 'ref_file')
    analysisflow.connect(inv_fnirt_xfm, 'inverse_warp', brain_warp, 'field_file')
    analysisflow.connect(inputspec, 'diff2highres',brain_warp,'postmat')

    analysisflow.connect(brain_warp, 'out_file',thresholdbin,'in_file')

    analysisflow.connect(thresholdbin,'out_file',outputspec, 'output_ventricle')

    return analysisflow

