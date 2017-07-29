function y = bandPassFilter15G35G(x)
%BANDPASSFILTER15G35G Filters input x and returns output y.

% MATLAB Code
% Generated by MATLAB(R) 9.1 and the DSP System Toolbox 9.3.
% Generated on: 29-Jul-2017 10:05:38

%#codegen

% To generate C/C++ code from this function use the codegen command. Type
% 'help codegen' for more information.

persistent Hd;

if isempty(Hd)
    
    % The following code was used to design the filter coefficients:
    % % Equiripple Bandpass filter designed using the FIRPM function.
    %
    % % All frequency values are in GHz.
    % Fs = 600;  % Sampling Frequency
    %
    % Fstop1 = 15;              % First Stopband Frequency
    % Fpass1 = 16;              % First Passband Frequency
    % Fpass2 = 34;              % Second Passband Frequency
    % Fstop2 = 35;              % Second Stopband Frequency
    % Dstop1 = 0.01;            % First Stopband Attenuation
    % Dpass  = 0.057501127785;  % Passband Ripple
    % Dstop2 = 0.01;            % Second Stopband Attenuation
    % dens   = 20;              % Density Factor
    %
    % % Calculate the order from the parameters using FIRPMORD.
    % [N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
    %                           0], [Dstop1 Dpass Dstop2]);
    %
    % % Calculate the coefficients using the FIRPM function.
    % b  = firpm(N, Fo, Ao, W, {dens});
    
    Hd = dsp.FIRFilter( ...
        'Numerator', [-0.00328662612537043 0.00304400411355789 ...
        0.00100889182085265 -0.000438850868031051 -0.00132696181743952 ...
        -0.00175190210005007 -0.00183121533693241 -0.00166979884539573 ...
        -0.00135063007552091 -0.000930213854910791 -0.000447793739480768 ...
        7.01502854637626e-05 0.000599082943233155 0.00111499796016897 ...
        0.00158879370188017 0.00199042106393222 0.00228805569576622 ...
        0.00245534154963468 0.00247186695365886 0.00232996454314494 ...
        0.00203343613493045 0.00160156054722416 0.00106464066323764 ...
        0.000464789604579942 -0.000151235860772574 -0.000733304467499819 ...
        -0.00123606394057572 -0.00162040071692578 -0.00185954669787317 ...
        -0.00194133868297838 -0.00186455281350489 -0.00165290903213645 ...
        -0.00132321370172661 -0.000928389644096508 -0.000498549297996215 ...
        -8.06686137619136e-05 0.000281683148107482 0.00056663832843261 ...
        0.000752974827932438 0.000831017438674074 0.000810340132873889 ...
        0.000712108483391345 0.000558954402278686 0.000381800229405877 ...
        0.000214905085461564 8.86016232063873e-05 2.15970925174225e-05 ...
        2.54452473931457e-05 0.000102114691495288 0.000243516189945656 ...
        0.000428414596537098 0.000629597193713792 0.000817087660524844 ...
        0.000963726196764916 0.00104440356130544 0.00104141933081865 ...
        0.000946025734548229 0.000762991498849653 0.00050681707540784 ...
        0.000201879978037315 -0.000124802121067749 -0.000439850320229998 ...
        -0.000715518199022344 -0.000921489111515442 -0.00104369214945667 ...
        -0.00107158841536745 -0.00101043550147316 -0.000874312333930782 ...
        -0.000682850390710184 -0.000463496237267545 -0.000245907062746899 ...
        -5.57654638128763e-05 8.33016004463481e-05 0.000157298059727724 ...
        0.000163417652851349 0.000108357570736643 6.22130345541688e-06 ...
        -0.000120004985850942 -0.000244338584990469 -0.000340856849165411 ...
        -0.000386727890652857 -0.00036413852768848 -0.000265911360153111 ...
        -9.56167729300418e-05 0.00013299819207584 0.000398542481372163 ...
        0.00067258959092237 0.000924216843429986 0.00112389144455689 ...
        0.00124726186220718 0.00127750450535923 0.00120701485317963 ...
        0.00104122818991195 0.000795143939641235 0.000494959911938399 ...
        0.0001698439572755 -0.000146837603193026 -0.00042401808134097 ...
        -0.00063674053344249 -0.000767259029494064 -0.000810149655534392 ...
        -0.000769857026063082 -0.000660550272112195 -0.00050601097637171 ...
        -0.000334925300521912 -0.000176626634467005 -5.82744976658197e-05 ...
        3.54460557013686e-07 -1.0857967414144e-05 -9.21309699669654e-05 ...
        -0.000232653641521895 -0.000411198658250568 -0.000600092633717949 ...
        -0.000768803766872346 -0.000887099008878146 -0.000929571375080824 ...
        -0.000879495148229897 -0.000730829839378788 -0.000489545120836011 ...
        -0.00017394567758044 0.000187299688457452 0.000559814519770569 ...
        0.000906359232992611 0.00119226879218583 0.00138778487208707 ...
        0.00147456426552143 0.00144529719748095 0.00130623335948495 ...
        0.00107575524605721 0.000781767039083333 0.000459576622214854 ...
        0.000145261220612687 -0.000126929976769509 -0.000329263874046824 ...
        -0.000444773108831045 -0.000469003181068826 -0.000410179786075459 ...
        -0.000288669487758325 -0.000132589537973004 2.48016638794034e-05 ...
        0.00014978004460669 0.000213283290565898 0.000195421764829461 ...
        8.71277846385453e-05 -0.000106999435333293 -0.000369435691282278 ...
        -0.000671492279966391 -0.000977180276253224 -0.00124677778594238 ...
        -0.00144314329747831 -0.00153534071431914 -0.00150355828780089 ...
        -0.0013411002439002 -0.00105769384191543 -0.000676102860823836 ...
        -0.000231686201338911 0.000233750105774306 0.00067460076808125 ...
        0.0010489608712587 0.00132247220935748 0.0014726379632844 ...
        0.0014918442985361 0.0013879139476711 0.00118267382896355 ...
        0.000909705636970838 0.000608861652961009 0.000322446166820045 ...
        8.83327106037182e-05 -6.46007754473701e-05 -0.000120022766340886 ...
        -7.60534235484649e-05 5.34057589958526e-05 0.00024153882750956 ...
        0.000451323732411049 0.000641182562493747 0.000770531430781153 ...
        0.00080563032052963 0.000723779343196941 0.000517819100199632 ...
        0.000196903400666209 -0.000213967619604631 -0.00067631617842782 ...
        -0.00114249206462973 -0.00156179402502924 -0.00188685197792869 ...
        -0.00207879994291101 -0.00211305718860965 -0.00198202892414378 ...
        -0.00169656602917932 -0.00128494133578885 -0.000789139792058816 ...
        -0.000260405419926899 0.000247950107186698 0.000685677229219765 ...
        0.00101256417908609 0.00120335255314622 0.00125050509943123 ...
        0.00116487981070605 0.000974234375788531 0.000719014140887445 ...
        0.000447341441194237 0.000207766146104695 4.30101507412732e-05 ...
        -1.57985268864564e-05 4.62087485043226e-05 0.000225128356668538 ...
        0.000498987482952428 0.000829723540282721 0.00116819794307105 ...
        0.00146073463385174 0.00165621568149936 0.00171261386921623 ...
        0.00160350177593761 0.00132151003027128 0.000880630533425065 ...
        0.000314727012247145 -0.000325599529563242 -0.000979666643054652 ...
        -0.00158276710096909 -0.00207436911658404 -0.0024056987236634 ...
        -0.00254541670913383 -0.00248380158302246 -0.00223355376173247 ...
        -0.00182877603090337 -0.00132022093085801 -0.000768922931848318 ...
        -0.000238756497498986 0.000211920534134381 0.000537653414595742 ...
        0.000711189884544474 0.000727381714018799 0.000603282292217724 ...
        0.000375881815823385 9.69427135301664e-05 -0.000173970560383575 ...
        -0.000377969785181649 -0.000464815686752683 -0.000400735933118763 ...
        -0.000172881036955175 0.000207907956725917 0.000707701528353655 ...
        0.00127274449972763 0.00183672844881298 0.00232804795130597 ...
        0.00267970825643838 0.00283733267352425 0.00276717019629157 ...
        0.00246010221506223 0.00193452313893483 0.00123401569727873 ...
        0.000422980845680688 -0.000420837060331192 -0.00121561355803282 ...
        -0.00188539489939478 -0.00236893007829313 -0.00262752023472654 ...
        -0.00264918621014181 -0.00245042682700779 -0.0020730550930257 ...
        -0.00157938864411663 -0.00104368684900553 -0.000542395231345026 ...
        -0.00014377255489376 0.00010062270506645 0.000163502771310431 ...
        4.51122137159589e-05 -0.000226125753595805 -0.000597548127718081 ...
        -0.000998977176520116 -0.0013524214177768 -0.00158226734961679 ...
        -0.00162623896286873 -0.00144404397705335 -0.00102405104923416 ...
        -0.000386106235334779 0.00041985905750699 0.00131908586472107 ...
        0.00222018841315067 0.0030263453705315 0.00364703081292596 ...
        0.00400936762030455 0.00406704471117621 0.00380737078325123 ...
        0.00325254949650954 0.00245808925083192 0.0015059340455678 ...
        0.000494947934959459 -0.000471296419986641 -0.00129694437769511 ...
        -0.00190596381952195 -0.00225124306662603 -0.00232048902275608 ...
        -0.0021369599407738 -0.00175667897039204 -0.00126035621781129 ...
        -0.000742602597978786 -0.000298931173632661 -1.28238124556578e-05 ...
        5.65116440890238e-05 -0.000117731960555206 -0.000525102533704015 ...
        -0.0011184379686099 -0.00181884112412 -0.00252569627192773 ...
        -0.00312915625539638 -0.00352504514379195 -0.00362833734856 ...
        -0.00338569425110275 -0.00278351631238203 -0.00185159591880646 ...
        -0.000660940004152005 0.000683277009458143 0.00205365700422301 ...
        0.00331549878296931 0.0043426336009192 0.00503314466871851 ...
        0.0053219765461137 0.00518877879110746 0.00466094023464689 ...
        0.00381035656364015 0.00274459412427396 0.00159339058362346 ...
        0.000492306101669976 -0.000434484172512492 -0.00109029134125441 ...
        -0.00141828233332273 -0.00140901642164944 -0.00110175637917725 ...
        -0.000579344758606377 4.29460819637243e-05 0.000632593561256163 ...
        0.00105834759673971 0.00120851275134326 0.00100729156301277 ...
        0.000426638903775299 -0.000507259059305278 -0.00171534306289911 ...
        -0.00307319027223928 -0.00442551398932953 -0.00560452615568813 ...
        -0.00645110588174172 -0.00683520257424018 -0.00667364766750201 ...
        -0.00594222974120174 -0.00468111446862241 -0.00299207638046581 ...
        -0.00102792029613376 0.00102477345274295 0.00296690277347673 ...
        0.00461121674661687 0.00580465323956328 0.00644758109654794 ...
        0.00650614650551551 0.00601692029734922 0.00508205608295472 ...
        0.00385676411739505 0.00252908705556652 0.0012954057329596 ...
        0.00033392772734608 -0.000219167097864196 -0.000289414929416928 ...
        0.000124403235381756 0.000947754414287736 0.00203837211039999 ...
        0.00320288915273744 0.00422173572020611 0.00487779948276268 ...
        0.00498650633740616 0.00442271708676905 0.00314149327143227 ...
        0.00118870819865408 -0.0012985662094508 -0.00410420725440432 ...
        -0.00695501346666282 -0.00955179445989282 -0.0116041478867957 ...
        -0.0128667213042378 -0.0131703326327861 -0.0124459822775461 ...
        -0.0107365596111877 -0.00819567088522479 -0.00507247479622976 ...
        -0.00168454725398517 0.00161955783172027 0.00450182489481855 ...
        0.00667734616387153 0.0079510099151643 0.00824476822354623 ...
        0.00761019367806528 0.00622600033387734 0.00437814465879763 ...
        0.00242582493169667 0.000755676146625718 -0.000269317815613622 ...
        -0.000361724454472118 0.000645783038896334 0.00276769679977573 ...
        0.00585055603469276 0.00957728730912253 0.0134933354578465 ...
        0.017050851570732 0.019668218610711 0.0207980516183748 ...
        0.0199968851080923 0.016987814201718 0.0117092270055593 ...
        0.00434240225380422 -0.00468565127212598 -0.0147261908649202 ...
        -0.024958236294348 -0.0344602329905804 -0.0422961986994441 ...
        -0.0476092008608115 -0.0497113373163834 -0.0481616612307747 ...
        -0.0428229932331479 -0.0338909278290166 -0.0218910978167643 ...
        -0.0076443301765924 0.00779828985197979 0.0232445587989121 ...
        0.0374658740388958 0.0493077648744371 0.0577955816888893 ...
        0.0622243741345219 0.0622243741345219 0.0577955816888893 ...
        0.0493077648744371 0.0374658740388958 0.0232445587989121 ...
        0.00779828985197979 -0.0076443301765924 -0.0218910978167643 ...
        -0.0338909278290166 -0.0428229932331479 -0.0481616612307747 ...
        -0.0497113373163834 -0.0476092008608115 -0.0422961986994441 ...
        -0.0344602329905804 -0.024958236294348 -0.0147261908649202 ...
        -0.00468565127212598 0.00434240225380422 0.0117092270055593 ...
        0.016987814201718 0.0199968851080923 0.0207980516183748 ...
        0.019668218610711 0.017050851570732 0.0134933354578465 ...
        0.00957728730912253 0.00585055603469276 0.00276769679977573 ...
        0.000645783038896334 -0.000361724454472118 -0.000269317815613622 ...
        0.000755676146625718 0.00242582493169667 0.00437814465879763 ...
        0.00622600033387734 0.00761019367806528 0.00824476822354623 ...
        0.0079510099151643 0.00667734616387153 0.00450182489481855 ...
        0.00161955783172027 -0.00168454725398517 -0.00507247479622976 ...
        -0.00819567088522479 -0.0107365596111877 -0.0124459822775461 ...
        -0.0131703326327861 -0.0128667213042378 -0.0116041478867957 ...
        -0.00955179445989282 -0.00695501346666282 -0.00410420725440432 ...
        -0.0012985662094508 0.00118870819865408 0.00314149327143227 ...
        0.00442271708676905 0.00498650633740616 0.00487779948276268 ...
        0.00422173572020611 0.00320288915273744 0.00203837211039999 ...
        0.000947754414287736 0.000124403235381756 -0.000289414929416928 ...
        -0.000219167097864196 0.00033392772734608 0.0012954057329596 ...
        0.00252908705556652 0.00385676411739505 0.00508205608295472 ...
        0.00601692029734922 0.00650614650551551 0.00644758109654794 ...
        0.00580465323956328 0.00461121674661687 0.00296690277347673 ...
        0.00102477345274295 -0.00102792029613376 -0.00299207638046581 ...
        -0.00468111446862241 -0.00594222974120174 -0.00667364766750201 ...
        -0.00683520257424018 -0.00645110588174172 -0.00560452615568813 ...
        -0.00442551398932953 -0.00307319027223928 -0.00171534306289911 ...
        -0.000507259059305278 0.000426638903775299 0.00100729156301277 ...
        0.00120851275134326 0.00105834759673971 0.000632593561256163 ...
        4.29460819637243e-05 -0.000579344758606377 -0.00110175637917725 ...
        -0.00140901642164944 -0.00141828233332273 -0.00109029134125441 ...
        -0.000434484172512492 0.000492306101669976 0.00159339058362346 ...
        0.00274459412427396 0.00381035656364015 0.00466094023464689 ...
        0.00518877879110746 0.0053219765461137 0.00503314466871851 ...
        0.0043426336009192 0.00331549878296931 0.00205365700422301 ...
        0.000683277009458143 -0.000660940004152005 -0.00185159591880646 ...
        -0.00278351631238203 -0.00338569425110275 -0.00362833734856 ...
        -0.00352504514379195 -0.00312915625539638 -0.00252569627192773 ...
        -0.00181884112412 -0.0011184379686099 -0.000525102533704015 ...
        -0.000117731960555206 5.65116440890238e-05 -1.28238124556578e-05 ...
        -0.000298931173632661 -0.000742602597978786 -0.00126035621781129 ...
        -0.00175667897039204 -0.0021369599407738 -0.00232048902275608 ...
        -0.00225124306662603 -0.00190596381952195 -0.00129694437769511 ...
        -0.000471296419986641 0.000494947934959459 0.0015059340455678 ...
        0.00245808925083192 0.00325254949650954 0.00380737078325123 ...
        0.00406704471117621 0.00400936762030455 0.00364703081292596 ...
        0.0030263453705315 0.00222018841315067 0.00131908586472107 ...
        0.00041985905750699 -0.000386106235334779 -0.00102405104923416 ...
        -0.00144404397705335 -0.00162623896286873 -0.00158226734961679 ...
        -0.0013524214177768 -0.000998977176520116 -0.000597548127718081 ...
        -0.000226125753595805 4.51122137159589e-05 0.000163502771310431 ...
        0.00010062270506645 -0.00014377255489376 -0.000542395231345026 ...
        -0.00104368684900553 -0.00157938864411663 -0.0020730550930257 ...
        -0.00245042682700779 -0.00264918621014181 -0.00262752023472654 ...
        -0.00236893007829313 -0.00188539489939478 -0.00121561355803282 ...
        -0.000420837060331192 0.000422980845680688 0.00123401569727873 ...
        0.00193452313893483 0.00246010221506223 0.00276717019629157 ...
        0.00283733267352425 0.00267970825643838 0.00232804795130597 ...
        0.00183672844881298 0.00127274449972763 0.000707701528353655 ...
        0.000207907956725917 -0.000172881036955175 -0.000400735933118763 ...
        -0.000464815686752683 -0.000377969785181649 -0.000173970560383575 ...
        9.69427135301664e-05 0.000375881815823385 0.000603282292217724 ...
        0.000727381714018799 0.000711189884544474 0.000537653414595742 ...
        0.000211920534134381 -0.000238756497498986 -0.000768922931848318 ...
        -0.00132022093085801 -0.00182877603090337 -0.00223355376173247 ...
        -0.00248380158302246 -0.00254541670913383 -0.0024056987236634 ...
        -0.00207436911658404 -0.00158276710096909 -0.000979666643054652 ...
        -0.000325599529563242 0.000314727012247145 0.000880630533425065 ...
        0.00132151003027128 0.00160350177593761 0.00171261386921623 ...
        0.00165621568149936 0.00146073463385174 0.00116819794307105 ...
        0.000829723540282721 0.000498987482952428 0.000225128356668538 ...
        4.62087485043226e-05 -1.57985268864564e-05 4.30101507412732e-05 ...
        0.000207766146104695 0.000447341441194237 0.000719014140887445 ...
        0.000974234375788531 0.00116487981070605 0.00125050509943123 ...
        0.00120335255314622 0.00101256417908609 0.000685677229219765 ...
        0.000247950107186698 -0.000260405419926899 -0.000789139792058816 ...
        -0.00128494133578885 -0.00169656602917932 -0.00198202892414378 ...
        -0.00211305718860965 -0.00207879994291101 -0.00188685197792869 ...
        -0.00156179402502924 -0.00114249206462973 -0.00067631617842782 ...
        -0.000213967619604631 0.000196903400666209 0.000517819100199632 ...
        0.000723779343196941 0.00080563032052963 0.000770531430781153 ...
        0.000641182562493747 0.000451323732411049 0.00024153882750956 ...
        5.34057589958526e-05 -7.60534235484649e-05 -0.000120022766340886 ...
        -6.46007754473701e-05 8.83327106037182e-05 0.000322446166820045 ...
        0.000608861652961009 0.000909705636970838 0.00118267382896355 ...
        0.0013879139476711 0.0014918442985361 0.0014726379632844 ...
        0.00132247220935748 0.0010489608712587 0.00067460076808125 ...
        0.000233750105774306 -0.000231686201338911 -0.000676102860823836 ...
        -0.00105769384191543 -0.0013411002439002 -0.00150355828780089 ...
        -0.00153534071431914 -0.00144314329747831 -0.00124677778594238 ...
        -0.000977180276253224 -0.000671492279966391 -0.000369435691282278 ...
        -0.000106999435333293 8.71277846385453e-05 0.000195421764829461 ...
        0.000213283290565898 0.00014978004460669 2.48016638794034e-05 ...
        -0.000132589537973004 -0.000288669487758325 -0.000410179786075459 ...
        -0.000469003181068826 -0.000444773108831045 -0.000329263874046824 ...
        -0.000126929976769509 0.000145261220612687 0.000459576622214854 ...
        0.000781767039083333 0.00107575524605721 0.00130623335948495 ...
        0.00144529719748095 0.00147456426552143 0.00138778487208707 ...
        0.00119226879218583 0.000906359232992611 0.000559814519770569 ...
        0.000187299688457452 -0.00017394567758044 -0.000489545120836011 ...
        -0.000730829839378788 -0.000879495148229897 -0.000929571375080824 ...
        -0.000887099008878146 -0.000768803766872346 -0.000600092633717949 ...
        -0.000411198658250568 -0.000232653641521895 -9.21309699669654e-05 ...
        -1.0857967414144e-05 3.54460557013686e-07 -5.82744976658197e-05 ...
        -0.000176626634467005 -0.000334925300521912 -0.00050601097637171 ...
        -0.000660550272112195 -0.000769857026063082 -0.000810149655534392 ...
        -0.000767259029494064 -0.00063674053344249 -0.00042401808134097 ...
        -0.000146837603193026 0.0001698439572755 0.000494959911938399 ...
        0.000795143939641235 0.00104122818991195 0.00120701485317963 ...
        0.00127750450535923 0.00124726186220718 0.00112389144455689 ...
        0.000924216843429986 0.00067258959092237 0.000398542481372163 ...
        0.00013299819207584 -9.56167729300418e-05 -0.000265911360153111 ...
        -0.00036413852768848 -0.000386727890652857 -0.000340856849165411 ...
        -0.000244338584990469 -0.000120004985850942 6.22130345541688e-06 ...
        0.000108357570736643 0.000163417652851349 0.000157298059727724 ...
        8.33016004463481e-05 -5.57654638128763e-05 -0.000245907062746899 ...
        -0.000463496237267545 -0.000682850390710184 -0.000874312333930782 ...
        -0.00101043550147316 -0.00107158841536745 -0.00104369214945667 ...
        -0.000921489111515442 -0.000715518199022344 -0.000439850320229998 ...
        -0.000124802121067749 0.000201879978037315 0.00050681707540784 ...
        0.000762991498849653 0.000946025734548229 0.00104141933081865 ...
        0.00104440356130544 0.000963726196764916 0.000817087660524844 ...
        0.000629597193713792 0.000428414596537098 0.000243516189945656 ...
        0.000102114691495288 2.54452473931457e-05 2.15970925174225e-05 ...
        8.86016232063873e-05 0.000214905085461564 0.000381800229405877 ...
        0.000558954402278686 0.000712108483391345 0.000810340132873889 ...
        0.000831017438674074 0.000752974827932438 0.00056663832843261 ...
        0.000281683148107482 -8.06686137619136e-05 -0.000498549297996215 ...
        -0.000928389644096508 -0.00132321370172661 -0.00165290903213645 ...
        -0.00186455281350489 -0.00194133868297838 -0.00185954669787317 ...
        -0.00162040071692578 -0.00123606394057572 -0.000733304467499819 ...
        -0.000151235860772574 0.000464789604579942 0.00106464066323764 ...
        0.00160156054722416 0.00203343613493045 0.00232996454314494 ...
        0.00247186695365886 0.00245534154963468 0.00228805569576622 ...
        0.00199042106393222 0.00158879370188017 0.00111499796016897 ...
        0.000599082943233155 7.01502854637626e-05 -0.000447793739480768 ...
        -0.000930213854910791 -0.00135063007552091 -0.00166979884539573 ...
        -0.00183121533693241 -0.00175190210005007 -0.00132696181743952 ...
        -0.000438850868031051 0.00100889182085265 0.00304400411355789 ...
        -0.00328662612537043]);
end

y = step(Hd,double(x));


% [EOF]
