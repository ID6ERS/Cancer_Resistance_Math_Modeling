% function for extracting data from CSV files.
function [ncr,ncg,cr,cg,hours,ncx0,cx0]=extractData()
m=csvread('Updated_payoff_tables.csv');
cols=[8,9,6,7,4,5];
noCis(1:61,1:6)=m(123:183,cols);
noCis(1:61,7:14)=m(1:61,2:9);
withCis(1:61,1:6)=m(184:244,cols);
withCis(1:61,7:14)=m(62:122,2:9);
ncr=noCis(:,1:2:end); cr=withCis(:,1:2:end);
ncg=noCis(:,2:2:end); cg=withCis(:,2:2:end);
ncx0=[ncr(1,1:7),ncg(1,1:7)];
cx0=[cr(1,1:7),cg(1,1:7)];
hours=m(1:61,1);
