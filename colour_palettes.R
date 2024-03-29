# colour_palettes set

primary_colour_setting <- "#E1BE6A"
relapse_colour_setting <- "#40B0A6"
gain_colour_setting <- "firebrick1"
amp_colour_setting <- "firebrick3"
loss_colour_setting <- "dodgerblue1"
del_colour_setting <- "dodgerblue4"
resistant_colour_setting <- "#F8766D"
sensitive_colour_setting <- "#00BFC4"
prior_lines_1_colour_setting <- "#ffffcc"
prior_lines_2_colour_setting <- "#c2e699"
prior_lines_3_colour_setting <- "#78c679"
prior_lines_4_colour_setting <- "#238443"
prior_lines_5_colour_setting <- "#006837"
stage_1_colour_setting <- "#fef0d9"
stage_2_colour_setting <- "#fdcc8a"
stage_3_colour_setting <- "#fc8d59"
stage_4_colour_setting <- "#d7301f"
brca_colour_setting <- "#E66100"
non_brca_colour_setting <- "#5D3A9B"

colour_palettes <- list(arx_rlps=c(arx=primary_colour_setting,
                                   rlps=relapse_colour_setting),
                        diagnosis_relapse=c(diagnosis=primary_colour_setting,
                                            relapse=relapse_colour_setting),
                        primary_relapse=c(primary=primary_colour_setting,
                                          relapse=relapse_colour_setting),
                        amplification_deletion=c(gain=gain_colour_setting,
                                                 amplification=amp_colour_setting,
                                                 loss=loss_colour_setting,
                                                 deletion=del_colour_setting),
                        amp_del=c(gain=gain_colour_setting,
                                  amplification=amp_colour_setting,
                                  loss=loss_colour_setting,
                                  deletion=del_colour_setting),
                        AMP_DEL=c(GAIN=gain_colour_setting,
                                   AMP=amp_colour_setting,
                                   LOSS=loss_colour_setting,
                                   DEL=del_colour_setting),
                        resist_sensitive=c(resistant=resistant_colour_setting,
                                           sensitive=sensitive_colour_setting),
                        prior_lines=c("1"=prior_lines_1_colour_setting,
                                      "2"=prior_lines_2_colour_setting,
                                      "3"=prior_lines_3_colour_setting,
                                      "4"=prior_lines_4_colour_setting,
                                      "5+"=prior_lines_5_colour_setting),
                        stage=c("1"=stage_1_colour_setting,
                                "2"=stage_2_colour_setting,
                                "3"=stage_3_colour_setting,
                                "4"=stage_4_colour_setting),
                        brca_carriers=c("BRCA-mutant"=brca_colour_setting,
                                        "wildtype"=non_brca_colour_setting))

rm(list = ls()[grep("*_colour_setting",ls())])
