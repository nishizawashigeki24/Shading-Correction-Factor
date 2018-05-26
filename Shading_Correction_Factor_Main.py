
# coding: utf-8

# # 日よけ効果係数算定ツール  Main
# 
# - 一次目標：日よけ効果係数算定ツールのpython上での再現
# 

# ## X. 日よけ効果係数算定ツール本体＋入出力
# 

# ### X.1 本体プログラム

# In[3]:


""" 日よけ効果係数計算プログラム本体 """
import Shading_Correction_Factor_Modules as SCFModule

def Calc_ShadingCorrectionFactor(Path00, FileName00, FileName01, ClimateZone, NDT, etaID, Azimuth, WSSize ):

    """ 引数(の例) """    
#    Path00 = "./SCFConfig01/"  # 設定ファイルのあるパス
#    FileName00 = "Zone.csv"    # 地点データのファイル
#    FileName01 = "IncidentAngleCharacteristics.csv"   # 窓ガラスの入射角特性ファイル
#    ClimateZone = 6            # 地域区分(1～8地域, 他、ユニークなID設定可)
#    NDT = 6                    # 1時間の分割数,ツールの標準は6分割
#    Azimuth = "南"             # 窓面の方位を16方位か、角度(-180°< Azimuth <= 180°)で入力
#    etaID = 1                  # 入射角特性のデータセットのID
#    WSSize = [1.1, 2.1, 0.9, 1.05, 1.07, 0.88, 0.85, 0.98, 2.05, 1.02, 0.96, 0.92, 1.01, 0.97, 0.24, 0.28, 0.21, 0.2]
    # WSSize = [X1, X2, X3, X1yp, X1ym, X3yp, X3ym, Y1, Y2, Y3, Y1xp, Y1xm, Y3xp, Y3xm, Zxp, Zxm, Zyp, Zym]
    # WSSize:窓および日よけの寸法一式
    
    """ リストの初期設定 """    
    Nh = [0 for i in range(8761)]
    SCF01 = [[[0 for j in range(NDT)] for i in range(8761)] for h in range(2)]
    
    """ \Zone.csv から地点データの読み込み(D.1) """
        # Zone:地域区分, City:都市, Latitude:緯度, Longitude:経度, SRFileName:日射量ファイル名
        # HStart:暖房開始日, HEnd:暖房終了日, CStart:冷房開始日, CEnd:冷房終了日
    [Zone, City, Latitude, Longitude, SRFileName, HStart, HEnd, CStart, CEnd]             = SCFModule.input_Point(ClimateZone, Path00, FileName00)

    """ 気象データ読み込み(D.2) """
        # SRHour[i][j]：1時間間隔データ
        #               i=0：1/1 0時 ～ 8760：12/31 24時
        #               j=0：月日時刻の5桁or6桁表記, j=1：法線面直達日射量[kcal/(m2h)]
        #               j=2：水平面天空日射量[kcal/(m2h)], j=3：暖房冷房判定タグ(暖房期:1, 冷房期:2, 非空調期:0)
    SRHour = SCFModule.input_SRData(Path00, SRFileName, HStart, HEnd, CStart, CEnd)
    
    """ 窓面の方位(A.8) """
        # Azwj:窓面の方位(-180°< Azwj <= 180°)
    Azwj = SCFModule.calc_Azwj(Azimuth)
    
    """ 天空日射の効果係数(C.3) """    
        # gammayp:天空日射の効果係数, 天空の形態係数の2倍
    gammayp = SCFModule.calc_gammayp(WSSize)

    """ 反射日射の効果係数(C.5) """    
        # gammaym:反射日射の効果係数, 地面の形態係数の2倍
    gammaym = SCFModule.calc_gammaym(WSSize)
    
    """ 窓ガラスの入射角特性読み込み(D.5) """  
        # etamax:直達日射に対する入射角特性最大値(入射角0)
        # etaisr:天空・反射日射に対する入射角特性(遮蔽なしの場合)
        # etakk:入射角特性の係数, (cosθ)^0 ～ (cosθ)^7 の各項の係数 
    [etaID0, etamax, etaisr, etakk] = SCFModule.input_IncidentAngleCharacteristics(etaID, Path00, FileName01)    
    
    """ 窓面積算定 """  
        # Awj:窓面積[m2]
    Awj = WSSize[1] * WSSize[8]    # Awj = X2 * Y2
    
    """ +++++++++++++++++++++++++++++++++++ 1時間のループ【1回目】→ Nhを通しで計算するためだけ ++++ """       
    for Hour00 in range(8761):
        """ 日時の計算(A.1)【1回目】 """
            # Hour00:1年間の通しの時刻, 1/1 0時がHour00=0, 12/31 24時がHour00=8760
            # NDay:1/1を"1", 12/31を"365"とする年頭からの通しの日数
            # NHour:時刻(0時～23時),12/31のみ24時あり
        [NDay,NHour] = SCFModule.calc_NDayNHour(Hour00)

        """ 正時±30分で太陽が地平線上にある時間刻み数のカウント(D.3) """
            # Nh[Hour00]:正時±30分で太陽が地平線上にある時間刻み数のカウント数
            # NDT:1時間の分割数,ツールの標準は6分割
        Nh[Hour00] = SCFModule.calc_Nh(Latitude, Longitude, NDay, NHour, NDT)
    
    """ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 1時間のループ【2回目】 ++++ """       
    for Hour00 in range(8760):      # ← 12/31 24時は回さないのでHour00=8759がループの最後
        """ 日時の計算(A.1)【2回目】 """
        [NDay,NHour] = SCFModule.calc_NDayNHour(Hour00)    

        """ 赤緯の計算(A.2) """
            # delta_d:赤緯[deg]
        deltad = SCFModule.calc_deltad(NDay)

        """ 均時差の計算(A.3) """
            # eed:均時差[hour]
        eed = SCFModule.calc_eed(NDay)
                
        """ 「月」の計算(D.2) """
        Month = SCFModule.calc_Month(SRHour.values[Hour00][0])
            
        """ +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ (1/NDT)時間のループ ++++ """ 
            # MM:1時間の内の(1/NDT)間隔の順番, 正時がMM=0, MM=0～NDT-1
        for MM in range(NDT):
            """ 時刻TTの計算(A.1) """
                # TT:時間分割MM毎の時刻[hour]
            TT = SCFModule.calc_TT(NHour, NDT, MM)

            """ 時角の計算(A.4) """
                # Tdt:時角[deg]
            Tdt = SCFModule.calc_Tdt(Longitude, eed, TT)

            """ 太陽高度とその正弦,余弦の計算(A.5,A.6) """
                # hsdt:太陽高度hsdt[deg], sinh:hsdtの正弦, cosh:hsdtの余弦
            sinh = SCFModule.calc_sinh(Latitude, deltad, Tdt)
            cosh = SCFModule.calc_cosh(sinh)
            hsdt = SCFModule.calc_hsdt(cosh, sinh)
            
            """ 太陽方位角の計算(A.7) """
                # Azsdt:太陽方位角[deg]
            Azsdt = SCFModule.calc_Azsdt(Latitude, deltad, Tdt, sinh, cosh)
            
            """ 窓面の法線ベクトルと太陽位置とのなす水平面上の角度の計算(A.9) """
                # Azwjdt:窓面の法線ベクトルと太陽位置とのなす水平面上の角度[deg]
            Azwjdt = SCFModule.calc_Azwjdt(Azwj, Azsdt)
            
            """ (1/NDT)分割MM番目の法線面直達日射量, 水平面天空日射量(D.4) """
                # Sddhm:(1/NDT)分割MM番目における法線面直達日射量[kcal/m2]
                # Ssdhm:(1/NDT)分割MM番目における水平面天空日射量[kcal/m2]
            Sddhm = SCFModule.calc_Sdhm(MM, NDT, sinh, SRHour.values[Hour00][1]
                                      , SRHour.values[Hour00+1][1], Nh[Hour00], Nh[Hour00+1])
            Ssdhm = SCFModule.calc_Sdhm(MM, NDT, sinh, SRHour.values[Hour00][2]
                                      , SRHour.values[Hour00+1][2], Nh[Hour00], Nh[Hour00+1]) 
            
            """ 直達日射の入射角(D.6) """
                # costheta:直達入射の窓面への入射角の余弦  
            costheta = SCFModule.calc_costheta(Azwjdt, cosh)
   
            """ 直達日射に対する窓ガラスの入射角特性(D.6) """
                # etajdt:日付NDay,時刻TTにおける入射角特性
            etajdt = SCFModule.calc_etajdt(costheta, etakk)
                        
            """ 直達日射が窓に射す部分の面積の計算(B.7) """
                # Ax:直達日射が窓に射す部分の面積[m2]
            Ax = SCFModule.calc_Ax(WSSize, Azwjdt, hsdt)

            """ 日よけ効果係数算定式の時刻TTにおける分子分母(E.2) """
                # SCF00[0]:(2)式の分子への加算分(窓面積をかけた値として)     
                # SCF00[1]:(2)式の分母への加算分(窓面積をかけた値として)
            [SCF01[0][Hour00][MM],SCF01[1][Hour00][MM]]                 = SCFModule.calc_SCF00(Sddhm, Ssdhm, etajdt, etaisr, costheta, sinh
                                       , Awj, Ax, gammayp, gammaym)
         
    SCF = SCFModule.Output_ShadingCorrectionFactor(SRHour, NDT, SCF01)
    
    return SCF

