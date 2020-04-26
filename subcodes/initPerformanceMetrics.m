function pfm_mtx = initPerformanceMetrics(general, varout)

ac_mtx = initPfmMtx(general, varout, 'ac');
nsdp_mtx = initPfmMtx(general, varout, 'nsdp');
tir_mtx = initPfmMtx(general, varout, 'tir');
stoi_mtx = initPfmMtx(general, varout, 'stoi');
psm_mtx =  initPfmMtx(general, varout, 'psm');
psmt_mtx =  initPfmMtx(general, varout, 'psmt');

pfm_mtx.ac_mtx = ac_mtx;
pfm_mtx.nsdp_mtx = nsdp_mtx;
pfm_mtx.tir_mtx = tir_mtx;
pfm_mtx.stoi_mtx = stoi_mtx;
pfm_mtx.psm_mtx = psm_mtx;
pfm_mtx.psmt_mtx = psmt_mtx;

end