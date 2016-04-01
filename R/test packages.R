vardisc <- c("branch","maritalstatus","sex","education","homestatus","profession","otherbusinesstype","havebankmandiriaccount","accountbank","asuransi","maincoverage","brand","phone","first2homezip","first2officezip","flagsame4zip","flagsame5zip","jointincome","flagotherloan","coverage","flagaccountbank","gendermarital","residencemarital","profaccount","genderwork","homezipdistance","flagbalance","dpdsr","dptenor","professionwork","brandasset","professioninc")
varcont <- c("ntf","umur","child","distance_from_tfs_to_cust","lamatinggal","lamabekerja","totalincome","totalotr","dsr","tenor","dp","averagebalance","LTV")

IV(data = dev,gbi = gbi,wt = wt,vardisc = vardisc,varcont = varcont)

library(Scorecard)
devtools::session_info()
