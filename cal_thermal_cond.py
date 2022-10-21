def calculate_thermal_cond_fromEMD(jx, jy, jz, temp, vol, dt, max_corr_time,savefig=False):
    """
    EMDの熱流束から熱伝導率の計算

    Parameters
    ----------
    jx,jy,jz: numpy.array
        熱流束J(t)，単位 eV/ps/A^2
    temp: float
        温度T．単位 K
    vol: float
        体積V．単位 A^3
    dt: int
        熱流束のサンプリングステップ．単位 fs
    max_corr_time: int
        相関関数を計算する最大時間．単位 fs

    Returns
    -------
    frequencies: numpy.ndarray
        np.array([np.arange(0, max_corr_time, 1)*dt/1000, kappa_list])
        相関時間[ps]とその時間までの熱伝導率
    """
    from scipy.signal import correlate
    import numpy as np
    import matplotlib.pyplot as plt
    from tqdm import tqdm
    
    # Σ J(0)・J(t)
    JJ = (correlate(jx, jx)+correlate(jy, jy)+correlate(jz, jz))/3
    
    # <J(0)・J(t)>
    # correlate mode="full"なので真ん中から計算し、各要素について適切に平均
    center = (len(JJ)-1)/2
    result = JJ[::-1][int(center):]
    result /= np.arange(1, len(result)+1, 1)[::-1]
    if savefig is True:
            plt.figure(figsize=(7, 6))
            plt.plot(np.arange(0, len(result), 1)*dt/1000, result)
            plt.xlim(-10, 1000)
            plt.xlabel("Correlation time[ps]", fontsize=14)
            plt.ylabel(r"HCACF $\left\langle J(0)J(t)\right\rangle$", fontsize=14)
            plt.savefig("HCACF.png",dpi=100)
        
    # V/kT^2 ∫<J(0)・J(t)>dt
    kappa_list = []
    max_corr_time = int(max_corr_time/dt)
    V = vol
    eV = 1.602176634e-19
    kb = 1.380649e-23
    T = temp
    scale = (V/kb/T/T)*1e-2*eV*eV/1e-12/1e-10

    for i in np.arange(0, max_corr_time, 1):
        kappa_list.append(np.sum(result[:i])*scale)
        
    if savefig is True:
        plt.figure(figsize=(7, 6))
        plt.plot(np.arange(0, max_corr_time, 1)*dt/1000, kappa_list)
        plt.xlabel("Correlation time[ps]", fontsize=14)
        plt.ylabel(r"$\kappa$ [W/m/K]", fontsize=14)
        plt.savefig("kappa.png", dpi=100)
    
    return np.array([np.arange(0, max_corr_time, 1)*dt/1000, kappa_list])
