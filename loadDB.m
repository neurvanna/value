function [subject, date] =  loadDB(sessionID)


db(1).subject = 'AL021';
db(1).date = '2019-06-05'; %k2, k3 good 1 2
db(2).subject = 'AL021';
db(2).date = '2019-06-06'; %k1,k2,k3, zo good 3 4 5 6
db(3).subject = 'AL021';
db(3).date = '2019-06-07';%k1,k2 good
db(4).subject = 'AL022';
db(4).date = '2019-06-19';%k1,k2,k3,zo good
db(5).subject = 'MW003';
db(5).date = '2019-08-11';%k1,k2,zo good
db(6).subject = 'MW003';
db(6).date = '2019-08-12';%k1,k2 good
db(7).subject = 'AL029';
db(7).date = '2019-10-23';%k3,zo ok
db(9).subject = 'AL026';
db(9).date = '2019-11-01';%k1, k2
db(10).subject = 'AL026';
db(10).date = '2019-11-02';%k1, k2 good
db(11).subject = 'AL026';
db(11).date = '2019-11-03';%k1, k2
db(13).subject = 'AL026';
db(13).date = '2019-11-05';%k1, k2
db(15).subject = 'AL026';
db(15).date = '2019-11-07';%k1, k2 ok
db(17).subject = 'AL026';
db(17).date = '2019-11-13';%k ok
db(18).subject = 'AL016';
db(18).date = '2019-07-11';%k1,k2
%db(8).subject = 'AL028';
%db(8).date = '2019-10-29';%k2


%db(12).subject = 'AL028';
%db(12).date = '2019-11-04';%k2

%db(14).subject = 'AL026';
%db(14).date = '2019-11-06';%k1, k2

%     db(16).subject = 'AL026';
%     db(16).date = '2019-11-12';%k2

subject = db(sessionID).subject;
date = db(sessionID).date;
