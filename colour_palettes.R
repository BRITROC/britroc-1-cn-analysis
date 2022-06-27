# colour_palettes set

primary_colour_setting <- "#E1BE6A"
relapse_colour_setting <- "#40B0A6"
gain_colour_setting <- "firebrick1"
amp_colour_setting <- "firebrick3"
loss_colour_setting <- "dodgerblue1"
del_colour_setting <- "dodgerblue4"
resistant_colour_setting <- "#F8766D"
sensitive_colour_setting <- "#00BFC4"
prior_lines_colour_1_setting <- "#ffffcc"
prior_lines_colour_2_setting <- "#c2e699"
prior_lines_colour_3_setting <- "#78c679"
prior_lines_colour_4_setting <- "#238443"
  
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
                        prior_lines=c("1"=prior_lines_colour_1_setting,
                                      "2"=prior_lines_colour_2_setting,
                                      "3"=prior_lines_colour_3_setting,
                                      "4"=prior_lines_colour_4_setting))

rm(list = ls()[grep("*_colour_setting",ls())])
