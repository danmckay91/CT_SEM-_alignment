
//open image

setBatchMode(true)

BSEI_PATH= "\\\\CSEG_2\\erc\\Mapping_01\\Raw_XCT_Data\\MAPPING_SCANS\\20151208_HUTCH_1050_SDK_R1_PCG_1[2015-12-10_13.31.26]\\20151208_HUTCH_1050_SDK_R1_PCG_1_COR\\alignment\\chemical\\PCG1-A_overlay_ChSEM_scale_by_0.482_flip_horrizontally.tif"
STACK_PATH = "\\\\CSEG_2\\erc\\Mapping_01\\Raw_XCT_Data\\MAPPING_SCANS\\20151208_HUTCH_1050_SDK_R1_PCG_1[2015-12-10_13.31.26]\\20151208_HUTCH_1050_SDK_R1_PCG_1_COR\\by_eye_4_cropped\\";
SAVE_PATH = "\\\\CSEG_2\\erc\\Mapping_01\\Raw_XCT_Data\\MAPPING_SCANS\\20151208_HUTCH_1050_SDK_R1_PCG_1[2015-12-10_13.31.26]\\20151208_HUTCH_1050_SDK_R1_PCG_1_COR\\alignment\\fiji_out_25_eye_four\\"+i+".txt"
open("BSEI_PATH");
rename("BSEI")


list = getFileList(STACK_PATH);
list = Array.sort(list)
for (i=470;i < 535; i++) { // was for (i=1;i < list.length; i++)
	print(list[i]);
	open(stackdir+list[i]);	
	rename("CT_slice");
	run("Extract SIFT Correspondences", "source_image=BSEI target_image=CT_slice initial_gaussian_blur=2 steps_per_scale_octave=5 minimum_image_size=50 maximum_image_size=1500 feature_descriptor_size=4 feature_descriptor_orientation_bins=8 closest/next_closest_ratio=0.95 filter maximal_alignment_error=25 minimal_inlier_ratio=0.02 minimal_number_of_inliers=3 expected_transformation=Translation");
	selectWindow("CT_slice");
	run("Measure");
	saveAs("Measurements", SAVE_PATH );
	close("CT_slice");
	run("Clear Results");
}




