%-----------------------------------------------------------------------
% Job saved on 18-Mar-2026 18:25:42 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM25 (25.01.02)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.tools.longit.pairwise.vols1 = {'/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/anat/02-3D_T1_3D_T1_20190527141834_2.nii,1'};
matlabbatch{1}.spm.tools.longit.pairwise.vols2 = {'/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/anat/02-3D_T1_3D_T1_20210708110821_2.nii,1'};
matlabbatch{1}.spm.tools.longit.pairwise.tdif = 1;
matlabbatch{1}.spm.tools.longit.pairwise.noise = NaN;
matlabbatch{1}.spm.tools.longit.pairwise.wparam = [0 0 100 25 100];
matlabbatch{1}.spm.tools.longit.pairwise.bparam = 1000000;
matlabbatch{1}.spm.tools.longit.pairwise.write_avg = 1;
matlabbatch{1}.spm.tools.longit.pairwise.write_jac = 0;
matlabbatch{1}.spm.tools.longit.pairwise.write_div = 1;
matlabbatch{1}.spm.tools.longit.pairwise.write_def = 0;
matlabbatch{2}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Pairwise Longitudinal Registration: Midpoint Average', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','avg', '()',{':'}));
matlabbatch{2}.spm.spatial.preproc.channel.biasreg = 0.0001;
matlabbatch{2}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{2}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{2}.spm.spatial.preproc.tissue(1).tpm = {'/Users/samira/Documents/MATLAB/spm/tpm/TPM.nii,1'};
matlabbatch{2}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{2}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(2).tpm = {'/Users/samira/Documents/MATLAB/spm/tpm/TPM.nii,2'};
matlabbatch{2}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{2}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(3).tpm = {'/Users/samira/Documents/MATLAB/spm/tpm/TPM.nii,3'};
matlabbatch{2}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{2}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(4).tpm = {'/Users/samira/Documents/MATLAB/spm/tpm/TPM.nii,4'};
matlabbatch{2}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{2}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(5).tpm = {'/Users/samira/Documents/MATLAB/spm/tpm/TPM.nii,5'};
matlabbatch{2}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{2}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(6).tpm = {'/Users/samira/Documents/MATLAB/spm/tpm/TPM.nii,6'};
matlabbatch{2}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{2}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{2}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{2}.spm.spatial.preproc.warp.reg = [0 0 0.1 0.01 0.04];
matlabbatch{2}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{2}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{2}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{2}.spm.spatial.preproc.warp.write = [0 1];
matlabbatch{2}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{2}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                              NaN NaN NaN];
matlabbatch{3}.spm.util.imcalc.input(1) = cfg_dep('Segment: INU corrected (1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
matlabbatch{3}.spm.util.imcalc.input(2) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{3}.spm.util.imcalc.input(3) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
matlabbatch{3}.spm.util.imcalc.input(4) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
matlabbatch{3}.spm.util.imcalc.output = 'skulStripped_biasCorrected_avgT1T2';
matlabbatch{3}.spm.util.imcalc.outdir = {'/Users/samira/Desktop/fMRI/sub-3002498'};
matlabbatch{3}.spm.util.imcalc.expression = '(i2 + i3 + i4) .* i1';
matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{3}.spm.util.imcalc.options.mask = 0;
matlabbatch{3}.spm.util.imcalc.options.interp = 0;
matlabbatch{3}.spm.util.imcalc.options.dtype = 16;
matlabbatch{4}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Image Calculator: ImCalc Computed Image: skulStripped_biasCorrected_avgT1T2', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(2) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(3) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(4) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
matlabbatch{4}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{4}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{4}.spm.spatial.normalise.write.woptions.interp = 7;
matlabbatch{4}.spm.spatial.normalise.write.woptions.prefix = 'w';
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).data.presubphasemag.phase = {'/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/fmap/15-fmri_fieldmap_gre_field_mapping_fmri_fieldmap_gre_field_mapping_20190527141834_15_e2_ph.nii,1'};
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).data.presubphasemag.magnitude = {'/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/fmap/14-fmri_fieldmap_gre_field_mapping_fmri_fieldmap_gre_field_mapping_20190527141834_14_e1.nii,1'};
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.et = [4.92 7.38];
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.maskbrain = 1;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.blipdir = -1;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.tert = 35.77;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.epifm = 0;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.ajm = 0;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.uflags.method = 'Mark3D';
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.uflags.fwhm = 10;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.uflags.pad = 0;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.uflags.ws = 1;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.mflags.template = {'/Users/samira/Documents/MATLAB/spm/matlabbatch/C:\Users\Samira\Documents\FMRI\spm25\spm\toolbox\FieldMap\T1.nii'};
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.mflags.fwhm = 5;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.mflags.nerode = 2;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.mflags.ndilate = 4;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.mflags.thresh = 0.5;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).defaults.defaultsval.mflags.reg = 0.02;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).session.epi = {'/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,1'};
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).matchvdm = 1;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).sessname = 'session';
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).writeunwarped = 1;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).anat = '';
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(1).matchanat = 1;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).data.presubphasemag.phase = {'/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/fmap/16-fmri_fieldmap_gre_field_mapping_fmri_fieldmap_gre_field_mapping_20210708110821_16_e2_ph.nii,1'};
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).data.presubphasemag.magnitude = {'/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/fmap/15-fmri_fieldmap_gre_field_mapping_fmri_fieldmap_gre_field_mapping_20210708110821_15_e1.nii,1'};
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.et = [4.92 7.38];
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.maskbrain = 1;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.blipdir = -1;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.tert = 35.77;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.epifm = 0;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.ajm = 0;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.uflags.method = 'Mark3D';
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.uflags.fwhm = 10;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.uflags.pad = 0;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.uflags.ws = 1;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.mflags.template = {'/Users/samira/Documents/MATLAB/spm/toolbox/FieldMap/T1.nii'};
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.mflags.fwhm = 5;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.mflags.nerode = 2;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.mflags.ndilate = 4;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.mflags.thresh = 0.5;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).defaults.defaultsval.mflags.reg = 0.02;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).session.epi = {'/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,1'};
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).matchvdm = 1;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).sessname = 'session';
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).writeunwarped = 1;
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).anat = '';
matlabbatch{5}.spm.tools.fieldmap.calculatevdm.subj(2).matchanat = 1;
%%
matlabbatch{6}.spm.spatial.realignunwarp.data(1).scans = {
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,1'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,2'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,3'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,4'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,5'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,6'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,7'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,8'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,9'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,10'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,11'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,12'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,13'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,14'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,15'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,16'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,17'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,18'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,19'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,20'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,21'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,22'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,23'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,24'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,25'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,26'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,27'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,28'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,29'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,30'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,31'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,32'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,33'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,34'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,35'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,36'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,37'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,38'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,39'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,40'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,41'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,42'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,43'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,44'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,45'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,46'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,47'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,48'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,49'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,50'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,51'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,52'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,53'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,54'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,55'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,56'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,57'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,58'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,59'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,60'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,61'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,62'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,63'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,64'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,65'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,66'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,67'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,68'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,69'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,70'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,71'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,72'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,73'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,74'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,75'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,76'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,77'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,78'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,79'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,80'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,81'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,82'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,83'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,84'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,85'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,86'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,87'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,88'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,89'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,90'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,91'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,92'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,93'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,94'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,95'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,96'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,97'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,98'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,99'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,100'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,101'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,102'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,103'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,104'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,105'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,106'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,107'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,108'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,109'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,110'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,111'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,112'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,113'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,114'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,115'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,116'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,117'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,118'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,119'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,120'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,121'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,122'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,123'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,124'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,125'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,126'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,127'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,128'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,129'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,130'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,131'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,132'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,133'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,134'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,135'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,136'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,137'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,138'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,139'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,140'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,141'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,142'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,143'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,144'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,145'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,146'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,147'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,148'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,149'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,150'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,151'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,152'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,153'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,154'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,155'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,156'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,157'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,158'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,159'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,160'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,161'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,162'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,163'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,164'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,165'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,166'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,167'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,168'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,169'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,170'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,171'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,172'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,173'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,174'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,175'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,176'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,177'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,178'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,179'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,180'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,181'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,182'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,183'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,184'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,185'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,186'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,187'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,188'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,189'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,190'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,191'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,192'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,193'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,194'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,195'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,196'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,197'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,198'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,199'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,200'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,201'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,202'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,203'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,204'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,205'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,206'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,207'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,208'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,209'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,210'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,211'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,212'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,213'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,214'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,215'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,216'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,217'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,218'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,219'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,220'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,221'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,222'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,223'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,224'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,225'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,226'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,227'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,228'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,229'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,230'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,231'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,232'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,233'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,234'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,235'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,236'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,237'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,238'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,239'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,240'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,241'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,242'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,243'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,244'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,245'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,246'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,247'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,248'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,249'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,250'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,251'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,252'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,253'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,254'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,255'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,256'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,257'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,258'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,259'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,260'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,261'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,262'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,263'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,264'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,265'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,266'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,267'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,268'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,269'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,270'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,271'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,272'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,273'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,274'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,275'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,276'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,277'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,278'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,279'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,280'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,281'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,282'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,283'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,284'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,285'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,286'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,287'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,288'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,289'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,290'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,291'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,292'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,293'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,294'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,295'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,296'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,297'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,298'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,299'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t04/func/12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.nii,300'
                                                          };
%%
matlabbatch{6}.spm.spatial.realignunwarp.data(1).pmscan(1) = cfg_dep('Calculate VDM: Voxel displacement map (Subj 1, Session 1)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','vdmfile', '{}',{1}));
%%
matlabbatch{6}.spm.spatial.realignunwarp.data(2).scans = {
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,1'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,2'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,3'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,4'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,5'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,6'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,7'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,8'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,9'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,10'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,11'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,12'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,13'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,14'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,15'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,16'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,17'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,18'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,19'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,20'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,21'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,22'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,23'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,24'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,25'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,26'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,27'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,28'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,29'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,30'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,31'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,32'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,33'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,34'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,35'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,36'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,37'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,38'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,39'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,40'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,41'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,42'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,43'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,44'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,45'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,46'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,47'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,48'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,49'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,50'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,51'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,52'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,53'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,54'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,55'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,56'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,57'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,58'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,59'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,60'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,61'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,62'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,63'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,64'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,65'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,66'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,67'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,68'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,69'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,70'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,71'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,72'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,73'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,74'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,75'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,76'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,77'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,78'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,79'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,80'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,81'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,82'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,83'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,84'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,85'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,86'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,87'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,88'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,89'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,90'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,91'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,92'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,93'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,94'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,95'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,96'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,97'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,98'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,99'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,100'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,101'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,102'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,103'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,104'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,105'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,106'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,107'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,108'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,109'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,110'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,111'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,112'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,113'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,114'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,115'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,116'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,117'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,118'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,119'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,120'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,121'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,122'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,123'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,124'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,125'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,126'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,127'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,128'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,129'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,130'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,131'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,132'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,133'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,134'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,135'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,136'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,137'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,138'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,139'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,140'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,141'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,142'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,143'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,144'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,145'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,146'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,147'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,148'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,149'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,150'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,151'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,152'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,153'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,154'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,155'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,156'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,157'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,158'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,159'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,160'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,161'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,162'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,163'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,164'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,165'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,166'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,167'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,168'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,169'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,170'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,171'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,172'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,173'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,174'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,175'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,176'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,177'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,178'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,179'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,180'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,181'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,182'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,183'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,184'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,185'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,186'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,187'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,188'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,189'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,190'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,191'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,192'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,193'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,194'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,195'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,196'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,197'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,198'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,199'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,200'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,201'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,202'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,203'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,204'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,205'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,206'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,207'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,208'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,209'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,210'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,211'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,212'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,213'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,214'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,215'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,216'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,217'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,218'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,219'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,220'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,221'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,222'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,223'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,224'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,225'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,226'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,227'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,228'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,229'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,230'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,231'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,232'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,233'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,234'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,235'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,236'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,237'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,238'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,239'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,240'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,241'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,242'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,243'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,244'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,245'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,246'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,247'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,248'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,249'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,250'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,251'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,252'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,253'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,254'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,255'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,256'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,257'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,258'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,259'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,260'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,261'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,262'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,263'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,264'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,265'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,266'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,267'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,268'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,269'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,270'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,271'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,272'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,273'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,274'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,275'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,276'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,277'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,278'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,279'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,280'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,281'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,282'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,283'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,284'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,285'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,286'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,287'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,288'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,289'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,290'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,291'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,292'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,293'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,294'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,295'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,296'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,297'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,298'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,299'
                                                          '/Users/samira/Desktop/fMRI/sub-3002498/ses-t06/func/13-fmri-MemoryTask_fmri-MemoryTask_20210708110821_13.nii,300'
                                                          };
%%
matlabbatch{6}.spm.spatial.realignunwarp.data(2).pmscan(1) = cfg_dep('Calculate VDM: Voxel displacement map (Subj 2, Session 1)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{2}, '.','vdmfile', '{}',{1}));
matlabbatch{6}.spm.spatial.realignunwarp.eoptions.quality = 0.95;
matlabbatch{6}.spm.spatial.realignunwarp.eoptions.sep = 1.5;
matlabbatch{6}.spm.spatial.realignunwarp.eoptions.fwhm = 1;
matlabbatch{6}.spm.spatial.realignunwarp.eoptions.rtm = 1;
matlabbatch{6}.spm.spatial.realignunwarp.eoptions.einterp = 7;
matlabbatch{6}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
matlabbatch{6}.spm.spatial.realignunwarp.eoptions.weight = '';
matlabbatch{6}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{6}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{6}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{6}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch{6}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch{6}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 2;
matlabbatch{6}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch{6}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch{6}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
matlabbatch{6}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{6}.spm.spatial.realignunwarp.uwroptions.jm = 1;
matlabbatch{6}.spm.spatial.realignunwarp.uwroptions.rinterp = 7;
matlabbatch{6}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
matlabbatch{6}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch{6}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
matlabbatch{7}.spm.temporal.st.scans{1}(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
matlabbatch{7}.spm.temporal.st.scans{2}(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 2)', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','uwrfiles'));
matlabbatch{7}.spm.temporal.st.nslices = 41;
matlabbatch{7}.spm.temporal.st.tr = 2.5;
matlabbatch{7}.spm.temporal.st.ta = 2.4390243902439;
matlabbatch{7}.spm.temporal.st.so = [2.42 2.36 2.3 2.2375 2.1775 2.1175 2.0575 1.9975 1.935 1.875 1.815 1.755 1.695 1.635 1.5725 1.5125 1.4525 1.3925 1.3325 1.27 1.21 1.15 1.09 1.03 0.9675 0.9075 0.8475 0.7875 0.7275 0.665 0.605 0.545 0.485 0.425 0.365 0.3025 0.2425 0.1825 0.1225 0.0625 0];
matlabbatch{7}.spm.temporal.st.refslice = 1.21;
matlabbatch{7}.spm.temporal.st.prefix = 'a';
matlabbatch{8}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Image Calculator: ImCalc Computed Image: skulStripped_biasCorrected_avgT1T2', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{8}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
matlabbatch{8}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{8}.spm.spatial.coreg.estimate.other(2) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 2)', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{2}, '.','files'));
matlabbatch{8}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{8}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{8}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{8}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
matlabbatch{9}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{9}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
matlabbatch{9}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{9}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{9}.spm.spatial.normalise.write.woptions.interp = 7;
matlabbatch{9}.spm.spatial.normalise.write.woptions.prefix = 'w';
matlabbatch{10}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{10}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{10}.spm.spatial.smooth.dtype = 0;
matlabbatch{10}.spm.spatial.smooth.im = 0;
matlabbatch{10}.spm.spatial.smooth.prefix = 's';
matlabbatch{11}.spm.util.imcalc.input(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{11}.spm.util.imcalc.output = 'brain_mask.nii';
matlabbatch{11}.spm.util.imcalc.outdir = {'/Users/samira/Desktop/fMRI/sub-3002498'};
matlabbatch{11}.spm.util.imcalc.expression = '(i2 + i3 + i4) > 0.2';
matlabbatch{11}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{11}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{11}.spm.util.imcalc.options.mask = 0;
matlabbatch{11}.spm.util.imcalc.options.interp = 0;
matlabbatch{11}.spm.util.imcalc.options.dtype = 2;
