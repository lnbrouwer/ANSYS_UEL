clearvars 
close all
fnum = 0;

cd Test_0A
delete *.jpg
cd ..

cd Test_0B
delete *.jpg
cd ..

cd Test_1A
delete *.jpg
cd ..

cd Test_2A
delete *.jpg
cd ..

cd Test_2B
delete *.jpg
cd ..

cd Test_2C
delete *.jpg
cd ..

cd Test_3A\Test_3Aa
delete *.jpg
cd ..
cd Test_3Ab
delete *.jpg
cd ..
cd Test_3Ac
delete *.jpg
cd .. 
cd ..


cd Test_3B\Test_3Ba
delete *.jpg
cd ..
cd Test_3Bb
delete *.jpg
cd ..
cd Test_3Bc
delete *.jpg
cd .. 
cd ..

cd Test_4A
delete *.jpg
cd ..

cd Test_4B
delete *.jpg
cd ..

cd Test_4C
delete *.jpg
cd ..

cd Test_4D
delete *.jpg
cd ..

cd Test_4E
delete *.jpg
cd ..


cd Test_5A
delete *.jpg
cd ..

cd Test_5B
delete *.jpg
cd ..

cd Test_5C
delete *.jpg
cd ..

cd Test_5D
delete *.jpg
cd ..


cd Test_6A
delete *.jpg
cd ..

cd Test_6B
delete *.jpg
cd ..




%test 0A
run .\Test_0A\compare_res.m
run .\Test_0A\compare_res_CC.m

%test 0B
run .\Test_0B\compare_res.m
run .\Test_0B\compare_res_CC.m

%test 1A
run .\Test_1A\compare_res.m
% run .\Test_1A\compare_res_CC.m

%test 2A
run .\Test_2A\compare_res.m
run .\Test_2A\compare_res_CC.m

%test 2B
run .\Test_2B\compare_res.m
run .\Test_2B\compare_res_CC.m

%test 2C
run .\Test_2C\compare_res.m
run .\Test_2C\compare_res_CC.m



%test 3A
run .\Test_3A\Test_3Aa\compare_res.m
run .\Test_3A\Test_3Aa\compare_res_CC.m
run .\Test_3A\Test_3Ab\compare_res.m
run .\Test_3A\Test_3Ab\compare_res_CC.m
run .\Test_3A\Test_3Ac\compare_res.m
run .\Test_3A\Test_3Ac\compare_res_CC.m


%test 3b
run .\Test_3B\Test_3Ba\compare_res.m
run .\Test_3B\Test_3Bb\compare_res.m
run .\Test_3B\Test_3Bc\compare_res.m
run .\Test_3B\Test_3Bc\compare_C_vs_b

%test 4A
run .\Test_4A\compare_res.m


%test 4B
run .\Test_4B\compare_res.m

%test 4C
run .\Test_4C\compare_res.m

%test 4D
run .\Test_4D\compare_res.m

%test 4E
run .\Test_4E\compare_res.m  %note that hotspot temp requires GUI
run .\Test_4E\compare_I_v_t.m



%test 5A
run .\Test_5A\compare_res.m
run .\Test_5A\compare_res_CC.m

%test 5B
run .\Test_5B\compare_res.m

%test 5C
run .\Test_5C\compare_res.m

%test 5D
run .\Test_5D\compare_res.m %note that hotspot temp requires GUI




%test 6A
run .\Test_6A\compare_res_cliq.m

%test 6B
run .\Test_6B\compare_res_cliq.m
run .\Test_6B\compare_res.m  %note that hotspot temp requires GUI

return 






