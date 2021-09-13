blocks <- list(one = c(307,278,279,277,276,305,306,296,297,298,299,269,270,151,150,152,153,149,145,146,147,148,229,230,224,223,225,226,211,216,215,214,213,212),
               two = c(120,121,119,124,123,122),
               three = c(125,142,143,141,144,140,139,138,2,1,137),
               four = c(107,108,131,4,109,111,110),
               five = c(160,161,164,163,162,165),
               six = c(316,315,314,159,158,157,156,90,89,87,86,84,85),
               seven = c(49,48,50,47,51,5,30,52,29,28,27,39,41,40,46,42,43,45,44,72,69,70,71),
               eight = c(292,293,294,302,272,301,273,274,300,271),
               nine = c(106,105,104,103,102,100,101,98,99),
               ten = c(26,25,23,24,22,63,76,57,21,55,56,77,53,54,78),
               eleven = c(235,232,233,234),
               twelve = c(79,82,81,80,83,154,155,130),
               thirteen = c(283,295,313,284,285,286,287,288),
               fourteen = c(309,310,311,312),
               fifteen = c(129,126,127,128),
               sixteen = c(308,280,97,96,95,94,92,91,93),
               seventeen = c(251,250,249,241,266,267,231),
               eighteen = c(14,13,204,201,202,228,227,203,218,12,11,59,210,195,217,196,6,7,9,8,20,18,19,16,58,17,15,60,10,61,197,62,198,199,200),
               nineteen = c(73,74,37,36,75,35,38,3,33,34,31,32,67,68,64,65,66,275,303,304),
               twenty = c(177,178,180,179,176),
               twenty_one = c(282,281,268,136,135,134,133,132,118,117),
               twenty_two = c(188,187,186,185,184,190,189,191,192,193),
               twenty_three = c(289,317,290,291),
               twenty_four = c(209,208,207,206,205,221,220,219,222,167,166,168,169,170,171,172,173,174,175),
               twenty_five = c(112,113,114,115,181,116,182,183,194),
               twenty_six = c(239,240,237,238,236,257,256,255,254,253,252),
               twenty_seven = c(242,243,244,245,246,247,248,263,262,264,265,258,259,260,261))



stud <- refine_data %>% dplyr::filter(Source_ID == "TW1_2006__Palomino")

for(i in 1:length(blocks)){
  stud <- stud %>% dplyr::mutate(Block = ifelse(Site_number %in% blocks[[i]], i, Block), SSB = ifelse(Site_number %in% blocks[[i]],paste("TW1_2006__Palomino 1",i), paste(SSB)))
}


refine_data <- refine_data %>% dplyr::filter(Source_ID != "TW1_2006__Palomino") %>% rbind(stud)



blocks <- list(one = c(65:68),
               two = c(94:100),
               three = c(71:93),
               four = c(51:59),
               five = c(60:64),
               six = c(69:70),
               seven = c(16:20,36:40),
               eight = c(3:4,23:24),
               nine = c(1:2,21:22),
               ten = c(11:15,31:35),
               eleven = c(41:45),
               twelve = c(46:50),
               thirteen = c(8:10,28:30),
               fourteen = c(5:7,25:27))


stud <- refine_data %>% dplyr::filter(Source_ID == "JB3_2019__Leso")

for(i in 1:length(blocks)){
  stud <- stud %>% dplyr::mutate(Block = ifelse(Site_number %in% blocks[[i]], i, Block), SSB = ifelse(Site_number %in% blocks[[i]],paste("JB3_2019__Leso 1",i), paste(SSB)))
}


refine_data <- refine_data %>% dplyr::filter(Source_ID != "JB3_2019__Leso") %>% rbind(stud)



stud <- refine_data %>% dplyr::filter(SS == "SE2_2013__Brandt 1")


stud$Block <- ifelse(stud$Site_number %in% c(67:69,5:7), "AAT", paste(stud$Block))


refine_data <- refine_data %>% dplyr::filter(SS != "SE2_2013__Brandt 1") %>% rbind(stud)




