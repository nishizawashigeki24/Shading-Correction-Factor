
# coding: utf-8

# # 日よけ効果係数算定ツール Modules
# 
# - 一次目標：日よけ効果係数算定ツールのpython上での再現
# 

# ## A. 窓まわりの設定(仕様書5)と太陽位置の算定(仕様書6.2)

# ### A.1 日時
# 
# - 日付`NDay`は1/1からの通しの日数
# - 時刻`NHour`は整数, 0時～23時が基本。12/31のみ24時あり
# - 時刻`TT`は、時間分割`MM`における時刻
#   - `MM`：1時間の内の$(1/$ `NDT`$)$間隔の順番, 正時が`MM`$=0$, `MM`$=0～$ `NDT`$-1$
# 
# 
# - 時刻`Hour01`は、前時刻の`MM` $=$ `NDT`$/2～$同時刻の`MM` $=$ `NDT`$/2-1$の時間分割を、その正時に属するものとして扱うための時刻 → 時間毎の効果係数を算定するために使用するだけ

# In[1]:


def calc_NDayNHour(Hour00):
    if Hour00==8760:
        [NDay, NHour] = [365, 24]
    else:
        NDay = Hour00//24 + 1
        NHour = Hour00 - (NDay - 1) * 24 

    return [NDay,NHour]


def calc_TT(NHour, NDT, MM):
    
    TT = NHour + MM / float(NDT)
    
    return TT


def calc_Hour01(TT):
    
    Hour01 = int(TT + 0.5)
    
    return Hour01


# ### A.2 赤緯の計算 (仕様書6.2 式(4))
# 
# - 赤緯$\delta_d [deg]$, $N$: 1月1日を$N=1$とした年頭からの通しの日数$[day]$
#   - 右辺の余弦のかっこ内の角度は$radian$単位となっているので注意
#   
# $$ \begin{align}
# \delta_d = (180 / \pi) & \{0.006322 - 0.405748 \cos (2 \pi N / 366 + 0.153231)\\
# & - 0.005880 \cos (4 \pi N / 366 + 0.207099)\\
# & - 0.003233 \cos (6 \pi N / 366 + 0.620129) \} \qquad \qquad \qquad (4) \\
# \end{align} $$

# In[2]:


""" 式(4) """
import numpy as np

def calc_deltad(NDay):
    deltad = (180 / np.pi) * (0.006322 - 0.405748 * np.cos(2 * np.pi * float(NDay) / 366 + 0.153231)
                                       - 0.005880 * np.cos(4 * np.pi * float(NDay) / 366 + 0.207099)
                                       - 0.003233 * np.cos(6 * np.pi * float(NDay) / 366 + 0.620129))
    
    return deltad


# ### A.3 均時差の計算 (仕様書6.2 式(6))
# 
# - 均時差$e_d[hour]$, $N$: 1月1日を$N=1$とした年頭からの通しの日数$[day]$
#   - 右辺の余弦のかっこ内の角度は$radian$単位となっているので注意
# 
# $$ \begin{align}
# e_d = -0.000279 &+ 0.122772 \cos (2 \pi N / 366 + 1.498311)\\
# & - 0.165458 \cos (4 \pi N / 366 - 1.261546)\\
# & - 0.005354 \cos (6 \pi N / 366 - 1.1571) \} \qquad \qquad \qquad (6) \\
# \end{align} $$

# In[3]:


""" 式(6) """
import numpy as np

def calc_eed(NDay):
    eed = ( -0.000279 + 0.122772 * np.cos(2 * np.pi * NDay / 366 + 1.498311)
                      - 0.165458 * np.cos(4 * np.pi * NDay / 366 - 1.261546)
                      - 0.005354 * np.cos(6 * np.pi * NDay / 366 - 1.1571)   )
    return eed


# ### A.4 時角の計算 (仕様書6.2 式(7))
# 
# - 時角$T_{d,t}[deg]$, 時刻$t[hour]$, 均時差$e_d[hour]$, 経度$L[deg]$
# 
# $$T_{d,t} = (t + e_d - 12) \times 15 + (L - 135)\qquad (7) $$

# In[4]:


""" 式(7) """
def calc_Tdt(Longitude, eed, TT):
    
    Tdt = ( TT + eed - 12) * 15 +(Longitude - 135)
    
    return Tdt


# ### A.5 太陽高度の正弦の計算 (仕様書6.2 式(8))
# 
# - 太陽高度$h_{S,d,t}[deg]$, 緯度$\phi[deg]$, 赤緯$\delta_d[deg]$,時角$T_{d,t}[deg]$
# 
# $$\sin h_{S,d,t} = max[0, \sin \phi \sin \delta_d + \cos \phi \cos \delta_d \cos T_{d,t}] \qquad (8) $$

# In[5]:


""" 式(8) """
import numpy as np

def calc_sinh(Latitude, deltad, Tdt):
    
    sinh = max(0, 
               np.sin(np.radians(Latitude)) * np.sin(np.radians(deltad)) 
                 + np.cos(np.radians(Latitude)) * np.cos(np.radians(deltad)) * np.cos(np.radians(Tdt)) )
   
    return sinh


# ### A.6 太陽高度とその余弦の計算 (仕様書6.2 式(9))
# 
# - 太陽高度$h_{S,d,t}[deg]$
# 
# $$\cos h_{S,d,t} = (1 - \sin ^2 h_{S,d,t})^{0.5} \qquad (9) $$
# 
# - 式(8), (9)より、$h_{S,d,t} = \tan^{-1} (\sin h_{S,d,t} / \cos h_{S,d,t})$

# In[6]:


""" 式(9)+α """
import numpy as np

def calc_cosh(sinh):
    
    cosh = (1 - sinh**2) **0.5
   
    return cosh


def calc_hsdt(cosh, sinh):
    
    hsdt = np.rad2deg(np.arctan( sinh / cosh ))
    
    return hsdt


# ### A.7 太陽方位角の計算 (仕様書6.2 式(10)～(12))
# 
# - 太陽方位角$A_{ZS,d,t}[deg]$, 太陽高度$h_{S,d,t}[deg]$, 赤緯$\delta_d[deg]$, 時角$T_{d,t}[deg]$, 緯度$\phi[deg]$
# 
# $$\sin A_{ZS,d,t} = \cos \delta_d \sin T_{d,t} / \cos h_{S,d,t} \qquad (10) $$
# 
# $$\cos A_{ZS,d,t} = (\sin h_{S,d,t} \sin \phi - \sin \delta_d) / (\cos h_{S,d,t} \cos \phi) \qquad (11) $$
# 
# $$
# A_{ZS,d,t} = \left\{
# \begin{array}{ll}
# \tan^{-1} (\sin A_{ZS,d,t} / \cos A_{ZS,d,t}) + 180 \hspace{24pt} (\sin A_{ZS,d,t} > 0, \cos A_{ZS,d,t} < 0)
# \\
# \tan^{-1} (\sin A_{ZS,d,t} / \cos A_{ZS,d,t}) - 180 \hspace{24pt} (\sin A_{ZS,d,t} < 0, \cos A_{ZS,d,t} < 0)
# \\
# 90 \hspace{136pt} (\sin A_{ZS,d,t} = 1, \cos A_{ZS,d,t} = 0)
# \\
# -90 \hspace{130pt} (\sin A_{ZS,d,t} = -1, \cos A_{ZS,d,t} = 0)
# \\
# \tan^{-1} (\sin A_{ZS,d,t} / \cos A_{ZS,d,t}) \hspace{48pt} (other)
# \end{array}
# \right.  \qquad (12) 
# $$
# 

# In[7]:


""" 式(10)～(12) """
import numpy as np

def calc_Azsdt(Latitude, deltad, Tdt, sinh, cosh):
    
    sinAzsdt = np.cos(np.radians(deltad)) * np.sin(np.radians(Tdt)) / cosh
    cosAzsdt = ( ( sinh * np.sin(np.radians(Latitude)) - np.sin(np.radians(deltad)) ) 
             /   ( cosh * np.cos(np.radians(Latitude)) ) )
    if abs(sinAzsdt) == 1:
        Azsdt = 90 * sinAzsdt
    elif sinAzsdt > 0 and cosAzsdt < 0:
        Azsdt = np.rad2deg(np.arctan( sinAzsdt / cosAzsdt )) + 180
    elif sinAzsdt < 0 and cosAzsdt < 0:        
        Azsdt = np.rad2deg(np.arctan( sinAzsdt / cosAzsdt) ) - 180   
    else:
        Azsdt = np.rad2deg(np.arctan( sinAzsdt / cosAzsdt) )
   
    return Azsdt


# ### A.8 窓面の方位 (仕様書5.2 図4)
# 
# - 窓面の方位は、以下の通り
#   - 北北東：$-157.5°$, 北東：$-135°$, …, 東：$-90°$, …, 南：$0°$, …, 西：$+90°$, …,北：$+180°$
#   - 角度指定も可：$-180°< A_{ZW,j} \leq +180°$
#   - デフォルトは8方位指定

# In[8]:


""" 窓面の方位 (仕様書5.2 図4) """
import sys

def calc_Azwj(Azimuth):
    
    Azimuth00 = ["北北東", "北東", "東北東", "東", "東南東", "南東", "南南東", "南"
                 , "南南西", "南西", "西南西", "西", "西北西", "北西", "北北西", "北" ]
    if Azimuth in Azimuth00:
        Azwj = (Azimuth00.index(Azimuth) - 7) * 22.5
    elif -180 < float(Azimuth) <= 180:
        Azwj = float(Azimuth) 
    else:
        sys.exit("窓面方位の入力が不適切です")        
    
    return Azwj


# ### A.9 窓面の法線ベクトルと太陽位置とのなす水平面上の角度の計算 (仕様書6.2 式(1))
# 
# - 窓面の法線ベクトルと太陽位置とのなす水平面上の角度$A_{ZW,j,d,t}[deg]$, 太陽方位角$A_{ZS,d,t}[deg]$, 外壁$j$の方位角$A_{ZW,j}[deg]$
# 
# $$
# A_{ZW,j,d,t} = \left\{
# \begin{array}{ll}
# A_{ZS,d,t} - A_{ZW,j} \hspace{48pt} (-180 < A_{ZS,d,t} - A_{ZW,j} \leq 180)
# \\
# A_{ZS,d,t} - A_{ZW,j} + 360 \hspace{24pt} (A_{ZS,d,t} - A_{ZW,j} \leq -180)
# \\
# A_{ZS,d,t} - A_{ZW,j} - 360 \hspace{24pt} (A_{ZS,d,t} - A_{ZW,j} \geq 180)
# \end{array}
# \right.  \qquad (1) 
# $$

# In[9]:


""" 式(1) """

def calc_Azwjdt(Azwj, Azsdt):
    
    Azwjdt = Azsdt - Azwj
    if Azwjdt < -180:
        Azwjdt += 360
    elif Azwjdt > 180:
        Azwjdt -= 360
        
    return Azwjdt


# ### A.10 窓まわり寸法のデータの持たせ方デフォルト (仕様書5.1 図2)

# In[10]:


import sys

def set_WSSize(WSSize1):
    
    WSSizeDict = {}
    
    for i in range(0, len(WSSize1), 2):
        WSSizeDict[WSSize1[i]] = WSSize1[i+1]

    # X1
    if "X1" not in WSSizeDict:
        WSSizeDict["X1"] = 0
    elif WSSizeDict["X1"] < 0:
        sys.exit("寸法X1の設定が不適切です")
    elif WSSizeDict["X1"] == "":
        WSSizeDict["X1"] = 0        
        
    # X2
    if "X2" not in WSSizeDict:
        sys.exit("寸法X2が設定されていません")
    elif WSSizeDict["X2"] <= 0 or WSSizeDict["X2"] == "":
        sys.exit("寸法X2の設定が不適切です")            
        
    # X3
    if "X3" not in WSSizeDict:
        WSSizeDict["X3"] = 0
    elif WSSizeDict["X3"] < 0:
        sys.exit("寸法X3の設定が不適切です")
    elif WSSizeDict["X3"] == "":
        WSSizeDict["X3"] = 0 
        
    # Y1
    if "Y1" not in WSSizeDict:
        WSSizeDict["Y1"] = 0
    elif WSSizeDict["Y1"] < 0:
        sys.exit("寸法Y1の設定が不適切です")
    elif WSSizeDict["Y1"] == "":
        WSSizeDict["Y1"] = 0 
        
    # Y2
    if "Y2" not in WSSizeDict:
        sys.exit("寸法Y2が設定されていません")
    elif WSSizeDict["Y2"] <= 0 or WSSizeDict["Y2"] == "":
        sys.exit("寸法Y2の設定が不適切です")            

    # Y3
    if "Y3" not in WSSizeDict:
        WSSizeDict["Y3"] = 0
    elif WSSizeDict["Y3"] < 0:
        sys.exit("寸法Y3の設定が不適切です")            
    elif WSSizeDict["Y3"] == "":
        WSSizeDict["Y3"] = 0 
        
    # Zxp
    if "Zxp" not in WSSizeDict:
        WSSizeDict["Zxp"] = 0
    elif WSSizeDict["Zxp"] < 0:
        sys.exit("寸法Zxpの設定が不適切です")
    elif WSSizeDict["Zxp"] == "":
        WSSizeDict["Zxp"] = 0 
        
    # Zxm
    if "Zxm" not in WSSizeDict:
        WSSizeDict["Zxm"] = 0
    elif WSSizeDict["Zxm"] < 0:
        sys.exit("寸法Zxmの設定が不適切です")           
    elif WSSizeDict["Zxm"] == "":
        WSSizeDict["Zxm"] = 0 
        
    # Zyp
    if "Zyp" not in WSSizeDict:
        WSSizeDict["Zyp"] = 0
    elif WSSizeDict["Zyp"] < 0:
        sys.exit("寸法Zypの設定が不適切です")     
    elif WSSizeDict["Zyp"] == "":
        WSSizeDict["Zyp"] = 0 
        
    # Zym
    if "Zym" not in WSSizeDict:
        WSSizeDict["Zym"] = 0
    elif WSSizeDict["Zym"] < 0:
        sys.exit("寸法Zymの設定が不適切です")
    elif WSSizeDict["Zym"] == "":
        WSSizeDict["Zym"] = 0 
        
    """ 以下、オプション扱い → 非入力は窓端から日よけの付け根までの距離とする """     
    # X1yp
    if "X1yp" not in WSSizeDict:
        WSSizeDict["X1yp"] = WSSizeDict["X1"] 
    elif WSSizeDict["X1yp"] < 0:
        sys.exit("寸法X1ypの設定が不適切です")        
    elif WSSizeDict["X1yp"] > WSSizeDict["X1"] or WSSizeDict["X1yp"] == "":
        WSSizeDict["X1yp"] = WSSizeDict["X1"]         
        
    # X1ym
    if "X1ym" not in WSSizeDict:
        WSSizeDict["X1ym"] = WSSizeDict["X1"] 
    elif WSSizeDict["X1ym"] < 0:
        sys.exit("寸法X1ymの設定が不適切です")
    elif WSSizeDict["X1ym"] > WSSizeDict["X1"] or WSSizeDict["X1ym"] == "":
        WSSizeDict["X1ym"] = WSSizeDict["X1"]       
        
    # X3yp
    if "X3yp" not in WSSizeDict:
        WSSizeDict["X3yp"] = WSSizeDict["X3"] 
    elif WSSizeDict["X3yp"] < 0:
        sys.exit("寸法X3ypの設定が不適切です")            
    elif WSSizeDict["X3yp"] > WSSizeDict["X3"] or WSSizeDict["X3yp"] == "":
        WSSizeDict["X3yp"] = WSSizeDict["X3"]   
        
    # X3ym
    if "X3ym" not in WSSizeDict:
        WSSizeDict["X3ym"] = WSSizeDict["X3"] 
    elif WSSizeDict["X3ym"] < 0:
        sys.exit("寸法X3ymの設定が不適切です")               
    elif WSSizeDict["X3ym"] > WSSizeDict["X3"] or WSSizeDict["X3ym"] == "":
        WSSizeDict["X3ym"] = WSSizeDict["X3"]   
        
    # Y1xp
    if "Y1xp" not in WSSizeDict:
        WSSizeDict["Y1xp"] = WSSizeDict["Y1"]
    elif WSSizeDict["Y1xp"] < 0:
        sys.exit("寸法Y1xpの設定が不適切です")     
    elif WSSizeDict["Y1xp"] > WSSizeDict["Y1"] or WSSizeDict["Y1xp"] == "":
        WSSizeDict["Y1xp"] = WSSizeDict["Y1"]  
        
    # Y1xm
    if "Y1xm" not in WSSizeDict:
        WSSizeDict["Y1xm"] = WSSizeDict["Y1"]  
    elif WSSizeDict["Y1xm"] < 0:
        sys.exit("寸法Y1xmの設定が不適切です")                
    elif WSSizeDict["Y1xm"] > WSSizeDict["Y1"] or WSSizeDict["Y1xm"] == "":
        WSSizeDict["Y1xm"] = WSSizeDict["Y1"]  
        
    # Y3xp
    if "Y3xp" not in WSSizeDict:
        WSSizeDict["Y3xp"] = WSSizeDict["Y3"]
    elif WSSizeDict["Y3xp"] < 0:
        sys.exit("寸法Y3xpの設定が不適切です")
    elif WSSizeDict["Y3xp"] > WSSizeDict["Y3"] or WSSizeDict["Y3xp"] == "":
        WSSizeDict["Y3xp"] = WSSizeDict["Y3"]
        
    # Y3xm
    if "Y3xm" not in WSSizeDict:
        WSSizeDict["Y3xm"] = WSSizeDict["Y3"]
    elif WSSizeDict["Y3xm"] < 0:
        sys.exit("寸法Y3xmの設定が不適切です")            
    elif WSSizeDict["Y3xm"] > WSSizeDict["Y3"] or WSSizeDict["Y3xm"] == "":
        WSSizeDict["Y3xm"] = WSSizeDict["Y3"]
        
    WSSize = [WSSizeDict["X1"],  WSSizeDict["X2"],  WSSizeDict["X3"],
              WSSizeDict["X1yp"],WSSizeDict["X1ym"],WSSizeDict["X3yp"],WSSizeDict["X3ym"],
              WSSizeDict["Y1"],  WSSizeDict["Y2"],  WSSizeDict["Y3"],
              WSSizeDict["Y1xp"],WSSizeDict["Y1xm"],WSSizeDict["Y3xp"],WSSizeDict["Y3xm"],
              WSSizeDict["Zxp"], WSSizeDict["Zxm"], WSSizeDict["Zyp"], WSSizeDict["Zym"] ]
            
    return WSSize


# ## B. 直達日射が窓に射す面積の計算 (仕様書6.3)

# ### B.1 太陽がx+側に位置する際のオーバーハングによる影の面積の算定式 (仕様書6.3.1 式(15))
# 
# $$
# A_{oh0+}(x,y) = \left\{
# \begin{array}{ll}
# 0 \hspace{24pt}(z_{y+}=0)
# \\
# \dfrac{1}{2} (x_{3y+} + x_2 / 2 - x) \dfrac{z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}}{z_{y+} \tan | A_{ZW,j,d,t} |} (x_{3y+} + x_2 / 2 - x)
# \\
# \hspace{30pt} \left( \begin{array}{ll}
# x_{3y+} + x_2 / 2 - x < z_{y+} \tan | A_{ZW,j,d,t} | \\
# y_{1} + y_2 / 2 - y \geq \dfrac{z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}}{z_{y+} \tan | A_{ZW,j,d,t} |} (x_{3y+} + x_2 / 2 - x) 
# \end{array} \right)
# \\
# \Bigl\{ (x_{3y+} + x_2 / 2 - x) - \dfrac{1}{2} (y_{1} + y_2 / 2 - y) \frac{z_{y+} \tan | A_{ZW,j,d,t} |}{z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}} \Bigr\} (y_{1} + y_2 / 2 - y)
# \\
# \hspace{30pt} \left( \begin{array}{ll}
# x_{3y+} + x_2 / 2 - x > \dfrac{z_{y+} \tan | A_{ZW,j,d,t} |}{z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}} (y_{1} + y_2 / 2 - y) \\
# y_{1} + y_2 / 2 - y < z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}
# \end{array} \right)
# \\
# ( x_{3y+} + x_2 / 2 - x - \dfrac{1}{2} z_{y+} \tan | A_{ZW,j,d,t} | ) z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}
# \\
# \hspace{30pt} \left( \begin{array}{ll}
# x_{3y+} + x_2 / 2 - x \geq z_{y+} \tan | A_{ZW,j,d,t} | \\
# y_{1} + y_2 / 2 - y \geq z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}
# \end{array} \right)
# \end{array}
# \right.  \qquad (15) 
# $$

# In[11]:


""" 式(15) """
import numpy as np

def calc_Aoh0p(XX, YY, WSSize, Azw, hs):
    [X1, X2, X3, X1yp, X1ym, X3yp, X3ym, Y1, Y2, Y3, Y1xp, Y1xm, Y3xp, Y3xm, Zxp, Zxm, Zyp, Zym] = WSSize
    
    X_th = X3yp + X2 / 2 - XX
    Y_th = Y1 + Y2 / 2 - YY
    X_th_Z = Zyp * np.tan(abs(np.radians(Azw)))  
    Y_th_Z = Zyp * np.tan(np.radians(hs)) / np.cos(np.radians(Azw))
        
    if X_th_Z == 0 or Y_th_Z <= 0:
        Aoh0p = 0    # 式(15)条件1 と 日よけが影を落とさない条件をあわせて処理
    else:
        Aoh0p = calc_Aoh0p00(X_th, Y_th, X_th_Z, Y_th_Z)
        
    return Aoh0p
        
    
def calc_Aoh0p00(X_th, Y_th, X_th_Z, Y_th_Z):
    
    if (X_th >= X_th_Z and Y_th >= Y_th_Z):
        Aoh0p00 = (X_th - X_th_Z / 2) * Y_th_Z       # 式(15)条件4
    elif Y_th * X_th_Z >= X_th * Y_th_Z:
        Aoh0p00 = X_th ** 2 * Y_th_Z / X_th_Z / 2       # 式(15)条件2
    else:
        Aoh0p00 = (X_th - Y_th / 2 * X_th_Z / Y_th_Z) * Y_th      # 式(15)条件3    
        
    return Aoh0p00


# ### B.2 太陽がx+側に位置する際のサイドフィンによる影の面積の算定式 (仕様書6.3.1 式(16))
# 
# $$
# A_{sf0+}(x,y) = \left\{
# \begin{array}{ll}
# 0 \hspace{24pt}(z_{x+}=0)
# \\
# \dfrac{1}{2} (y_{1x+} + y_2 / 2 - y) \dfrac{z_{x+} \tan | A_{ZW,j,d,t} |}{z_{x+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}} (y_{1x+} + y_2 / 2 - y)
# \\
# \hspace{30pt} \left( \begin{array}{ll}
# y_{1x+} + y_2 / 2 - y < z_{x+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t} \\
# x_{3} + x_2 / 2 - x \geq \dfrac{z_{x+} \tan | A_{ZW,j,d,t} |}{z_{x+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}} (y_{1x+} + y_2 / 2 - y) 
# \end{array} \right)
# \\
# \Bigl\{ (y_{1x+} + y_2 / 2 - y) - \dfrac{1}{2} (x_{3} + x_2 / 2 - x) \dfrac{z_{x+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}}{z_{x+} \tan | A_{ZW,j,d,t} |} \Bigr\} (x_{3} + x_2 / 2 - x)
# \\
# \hspace{30pt} \left( \begin{array}{ll}
# y_{1x+} + y_2 / 2 - y > \dfrac{z_{x+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}}{z_{x+} \tan | A_{ZW,j,d,t} |} (x_{3} + x_2 / 2 - x) \\
# x_{3} + x_2 / 2 - x < z_{x+} \tan | A_{ZW,j,d,t} |
# \end{array} \right)
# \\
# ( y_{1x+} + y_2 / 2 - y - \dfrac{1}{2} z_{x+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t} ) \; z_{x+} \tan | A_{ZW,j,d,t} |
# \\
# \hspace{30pt} \left( \begin{array}{ll}
# y_{1x+} + y_2 / 2 - y \geq z_{x+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t} \\
# x_{3} + x_2 / 2 - x \geq z_{x+} \tan | A_{ZW,j,d,t} |
# \end{array} \right)
# \end{array}
# \right.  \qquad (16) 
# $$
# 
# - コード中では、座標を入れ替えて、`calc_Aoh0p00`を叩くことで対応
#   - 式$(15)$の変数 → 式$(16)$の変数
#   - $x$ → $y$
#   - $x_2$ → $y_2$
#   - $x_{3y+}$ → $y_{1x+}$
#   - $y$ → $x$
#   - $y_1$ → $x_3$
#   - $y_2$ → $x_2$
#   - $z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}$ → $z_{x+} \tan | A_{ZW,j,d,t} |$
#   - $z_{y+} \tan | A_{ZW,j,d,t} |$ → $z_{x+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}$

# In[12]:


""" 式(16) """
import numpy as np

def calc_Asf0p(XX, YY, WSSize, Azw, hs):
    
    [X1, X2, X3, X1yp, X1ym, X3yp, X3ym, Y1, Y2, Y3, Y1xp, Y1xm, Y3xp, Y3xm, Zxp, Zxm, Zyp, Zym] = WSSize
    
    X_th = Y1xp + Y2 / 2 - YY
    Y_th = X3 + X2 / 2 - XX
    X_th_Z = Zxp * np.tan(np.radians(hs)) / np.cos(np.radians(Azw))  
    Y_th_Z = Zxp * np.tan(abs(np.radians(Azw)))
       
    if X_th_Z == 0 or Y_th_Z <= 0:
        Aoh0p = 0    # 式(16)条件1 と 日よけが影を落とさない条件をあわせて処理
    else:
        Aoh0p = calc_Aoh0p00(X_th, Y_th, X_th_Z, Y_th_Z)
        
    return Aoh0p


# ### B.3 太陽がx+側に位置する際の日射が射す部分の面積の計算式 (仕様書6.3.1 式(14))
# 
# $$ \begin{align}
# A_{wind,j,x+,d,t} &= (x_2 + x_3)(y_1 + y_2) - A_{oh0+}(-x_2 / 2, -y_2 / 2) - A_{sf0+}(-x_2 / 2, -y_2 / 2) \\
# &- \{ (x_2 + x_3) y_1 - A_{oh0+}(-x_2 / 2, y_2 / 2) - A_{sf0+}(-x_2 / 2, y_2 / 2) \} \\
# &- \{ x_3 (y_1 + y_2) - A_{oh0+}( x_2 / 2, -y_2 / 2) - A_{sf0+}( x_2 / 2, -y_2 / 2) \} \\
# &+ x_3 y_1 - A_{oh0+}( x_2 / 2, y_2 / 2) - A_{sf0+}( x_2 / 2, y_2 / 2) \qquad \qquad \qquad (14) \\
# \end{align} $$

# In[13]:


""" 式(14) """

def calc_Axp(WSSize, Azw, hs):
    
    [X1, X2, X3, X1yp, X1ym, X3yp, X3ym, Y1, Y2, Y3, Y1xp, Y1xm, Y3xp, Y3xm, Zxp, Zxm, Zyp, Zym] = WSSize
    
    if hs > 0 and -90 < Azw < 0:
        Axp = ( (X2 + X3) * (Y1 + Y2) 
               - calc_Aoh0p(-X2/2, -Y2/2, WSSize, Azw, hs) 
               - calc_Asf0p(-X2/2, -Y2/2, WSSize, Azw, hs) )\
            - ( (X2 + X3) * Y1        
               - calc_Aoh0p(-X2/2,  Y2/2, WSSize, Azw, hs) 
               - calc_Asf0p(-X2/2,  Y2/2, WSSize, Azw, hs) ) \
            - ( X3 * (Y1 + Y2)        
               - calc_Aoh0p( X2/2, -Y2/2, WSSize, Azw, hs) 
               - calc_Asf0p( X2/2, -Y2/2, WSSize, Azw, hs) ) \
            + ( X3 * Y1               
               - calc_Aoh0p( X2/2,  Y2/2, WSSize, Azw, hs) 
               - calc_Asf0p( X2/2,  Y2/2, WSSize, Azw, hs) )       
        Axp = max(0, min(Axp, X2 * Y2))    #負値は0に、X2*Y2を超える場合はX2*Y2で頭打ち
    else:
        Axp = 0
        
    return Axp


# ### B.4 太陽がx-側に位置する際のオーバーハングによる影の面積の算定式 (仕様書6.3.2 式(19))
# 
# $$
# A_{oh0-}(x,y) = \left\{
# \begin{array}{ll}
# 0 \hspace{24pt}(z_{y+}=0)
# \\
# \dfrac{1}{2} (x_{1y+} + x_2 / 2 + x) \dfrac{z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}}{z_{y+} \tan A_{ZW,j,d,t}} (x_{1y+} + x_2 / 2 + x)
# \\
# \hspace{30pt} \left( \begin{array}{ll}
# x_{1y+} + x_2 / 2 + x < z_{y+} \tan A_{ZW,j,d,t} \\
# y_{1} + y_2 / 2 - y \geq \dfrac{z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}}{z_{y+} \tan A_{ZW,j,d,t}} (x_{1y+} + x_2 / 2 + x) 
# \end{array} \right)
# \\
# \Bigl\{ (x_{1y+} + x_2 / 2 + x) - \dfrac{1}{2} (y_{1} + y_2 / 2 - y) \dfrac{z_{y+} \tan A_{ZW,j,d,t}}{z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}} \Bigr\} (y_{1} + y_2 / 2 - y)
# \\
# \hspace{30pt} \left( \begin{array}{ll}
# x_{1y+} + x_2 / 2 + x > \dfrac{z_{y+} \tan A_{ZW,j,d,t}}{z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}} (y_{1} + y_2 / 2 - y) \\
# y_{1} + y_2 / 2 - y < z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}
# \end{array} \right)
# \\
# ( x_{1y+} + x_2 / 2 + x - \dfrac{1}{2} z_{y+} \tan A_{ZW,j,d,t} ) z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}
# \\
# \hspace{30pt} \left( \begin{array}{ll}
# x_{1y+} + x_2 / 2 + x \geq z_{y+} \tan A_{ZW,j,d,t} \\
# y_{1} + y_2 / 2 - y \geq z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}
# \end{array} \right)
# \end{array}
# \right.  \qquad (19) 
# $$
# 
# - コード中では、座標を入れ替えて、`calc_Aoh0p00`を叩くことで対応
#   - 式$(15)$の変数 → 式$(19)$の変数
#   - $x$ → $-x$
#   - $x_{3y+}$ → $x_{1y+}$
#   - $A_{ZW,j,d,t}$ → $-A_{ZW,j,d,t}$

# In[14]:


""" 式(19) """
import numpy as np

def calc_Aoh0m(XX, YY, WSSize, Azw, hs):
    
    [X1, X2, X3, X1yp, X1ym, X3yp, X3ym, Y1, Y2, Y3, Y1xp, Y1xm, Y3xp, Y3xm, Zxp, Zxm, Zyp, Zym] = WSSize
    
    X_th = X1yp + X2 / 2 + XX
    Y_th = Y1 + Y2 / 2 - YY
    X_th_Z = Zyp * np.tan(abs(np.radians(Azw)))  
    Y_th_Z = Zyp * np.tan(np.radians(hs)) / np.cos(np.radians(Azw))
        
    if X_th_Z == 0 or Y_th_Z <= 0:
        Aoh0m = 0    # 式(19)条件1 と 日よけが影を落とさない条件をあわせて処理
    else:
        Aoh0m = calc_Aoh0p00(X_th, Y_th, X_th_Z, Y_th_Z)
        
    return Aoh0m


# ### B.5 太陽がx-側に位置する際のサイドフィンによる影の面積の算定式 (仕様書6.3.2 式(20))
# 
# $$
# A_{sf0-}(x, y) = \left\{
# \begin{array}{ll}
# 0 \hspace{24pt}(z_{x-}=0)
# \\
# \dfrac{1}{2} (y_{1x-} + y_2 / 2 - y) \dfrac{z_{x-} \tan A_{ZW,j,d,t}}{z_{x-} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}} (y_{1x-} + y_2 / 2 - y)
# \\
# \hspace{30pt} \left( \begin{array}{ll}
# y_{1x-} + y_2 / 2 - y < z_{x-} \tan h_{S,d,t} / \cos A_{ZW,j,d,t} \\
# x_{1} + x_2 / 2 + x \geq \dfrac{z_{x-} \tan A_{ZW,j,d,t}}{z_{x-} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}} (y_{1x-} + y_2 / 2 - y) 
# \end{array} \right)
# \\
# \Bigl\{ (y_{1x-} + y_2 / 2 - y) - \dfrac{1}{2} (x_{1} + x_2 / 2 + x) \dfrac{z_{x-} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}}{z_{x-} \tan A_{ZW,j,d,t} } \Bigr\} (x_{1} + x_2 / 2 + x)
# \\
# \hspace{30pt} \left( \begin{array}{ll}
# y_{1x-} + y_2 / 2 - y > \dfrac{z_{x-} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}}{z_{x-} \tan A_{ZW,j,d,t}} (x_{1} + x_2 / 2 + x) \\
# x_{1} + x_2 / 2 + x < z_{x-} \tan A_{ZW,j,d,t}
# \end{array} \right)
# \\
# ( y_{1x-} + y_2 / 2 - y - \dfrac{1}{2} z_{x-} \tan h_{S,d,t} / \cos A_{ZW,j,d,t} ) \; z_{x-} \tan A_{ZW,j,d,t}
# \\
# \hspace{30pt} \left( \begin{array}{ll}
# y_{1x-} + y_2 / 2 - y \geq z_{x-} \tan h_{S,d,t} / \cos A_{ZW,j,d,t} \\
# x_{1} + x_2 / 2 + x \geq z_{x-} \tan A_{ZW,j,d,t}
# \end{array} \right)
# \end{array}
# \right.  \qquad (20) 
# $$
# 
# - コード中では、座標を入れ替えて、`calc_Aoh0p00`を叩くことで対応
#   - 式$(15)$の変数 → 式$(16)$の変数 → 式$(20)$の変数
#   - $x$            → $y$            → $y$
#   - $x_2$          → $y_2$          → $y_2$
#   - $x_{3y+}$      → $y_{1x+}$      → $y_{1x-}$
#   - $y$            → $x$            → $-x$
#   - $y_1$          → $x_3$          → $x_1$
#   - $y_2$          → $x_2$          → $x_2$
#   - $z_{y+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}$ → $z_{x+} \tan | A_{ZW,j,d,t} |$ → $z_{x-} \tan A_{ZW,j,d,t}$
#   - $z_{y+} \tan | A_{ZW,j,d,t} |$ → $z_{x+} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}$ → $z_{x-} \tan h_{S,d,t} / \cos A_{ZW,j,d,t}$

# In[15]:


""" 式(20) """
import numpy as np

def calc_Asf0m(XX, YY, WSSize, Azw, hs):
    
    [X1, X2, X3, X1yp, X1ym, X3yp, X3ym, Y1, Y2, Y3, Y1xp, Y1xm, Y3xp, Y3xm, Zxp, Zxm, Zyp, Zym] = WSSize
    
    X_th = Y1xm + Y2 / 2 - YY
    Y_th = X1 + X2 / 2 + XX
    X_th_Z = Zxm * np.tan(np.radians(hs)) / np.cos(np.radians(Azw))  
    Y_th_Z = Zxm * np.tan(abs(np.radians(Azw)))
       
    if X_th_Z == 0 or Y_th_Z <= 0:
        Asf0m = 0    # 式(20)条件1 と 日よけが影を落とさない条件をあわせて処理
    else:
        Asf0m = calc_Aoh0p00(X_th, Y_th, X_th_Z, Y_th_Z)
        
    return Asf0m


# ### B.6 太陽がx-側に位置する際の日射が射す部分の面積の計算式 (仕様書6.3.2 式(18))
# 
# $$ \begin{align}
# A_{wind,j,x-,d,t} &= (x_1 + x_2)(y_1 + y_2) - A_{oh0-}( x_2 / 2, -y_2 / 2) - A_{sf0-}( x_2 / 2, -y_2 / 2) \\
# &- \{ (x_1 + x_2) y_1 - A_{oh0-}( x_2 / 2, y_2 / 2) - A_{sf0-}( x_2 / 2, y_2 / 2) \} \\
# &- \{ x_1 (y_1 + y_2) - A_{oh0-}(-x_2 / 2, -y_2 / 2) - A_{sf0-}(-x_2 / 2, -y_2 / 2) \} \\
# &+ x_1 y_1 - A_{oh0-}(-x_2 / 2, y_2 / 2) - A_{sf0-}(-x_2 / 2, y_2 / 2) \qquad \qquad \qquad (18) \\
# \end{align} $$

# In[16]:


""" 式(18) """

def calc_Axm(WSSize, Azw, hs):
    
    [X1, X2, X3, X1yp, X1ym, X3yp, X3ym, Y1, Y2, Y3, Y1xp, Y1xm, Y3xp, Y3xm, Zxp, Zxm, Zyp, Zym] = WSSize

    if hs > 0 and 0 <= Azw < 90:
        Axm = ( (X1 + X2) * (Y1 + Y2) 
               - calc_Aoh0m( X2/2, -Y2/2, WSSize, Azw, hs) 
               - calc_Asf0m( X2/2, -Y2/2, WSSize, Azw, hs) )\
            - ( (X1 + X2) * Y1        
               - calc_Aoh0m( X2/2,  Y2/2, WSSize, Azw, hs) 
               - calc_Asf0m( X2/2,  Y2/2, WSSize, Azw, hs) ) \
            - ( X1 * (Y1 + Y2)        
               - calc_Aoh0m(-X2/2, -Y2/2, WSSize, Azw, hs) 
               - calc_Asf0m(-X2/2, -Y2/2, WSSize, Azw, hs) ) \
            + ( X1 * Y1               
               - calc_Aoh0m(-X2/2,  Y2/2, WSSize, Azw, hs) 
               - calc_Asf0m(-X2/2,  Y2/2, WSSize, Azw, hs) )       
        Axm = max(0, min(Axm, X2 * Y2))    #負値は0に、X2*Y2を超える場合はX2*Y2で頭打ち
    else:
        Axm = 0
        
    return Axm


# ### B.7 直達日射が窓に射す部分の面積の計算
# 
# - 太陽高度$h_{S,d,t}[deg]$が、$h_{S,d,t}>0°$で計算
# - 窓面の法線ベクトルと太陽位置とのなす水平面上の角度$A_{ZW,j,d,t}$に応じて以下。
#   - $-90°<A_{ZW,j,d,t} < 0°$：$A_{wind,j,x-,d,t}$を計算($A_{wind,j,x+,d,t} = 0$)
#   - $0°\leq A_{ZW,j,d,t} < 90°$：$A_{wind,j,x+,d,t}$を計算($A_{wind,j,x-,d,t} = 0$)
# 

# In[17]:


def calc_Ax(WSSize, Azw, hs):

    if hs > 0 and -90 < Azw < 0:
        Ax = calc_Axp(WSSize, Azw, hs)
    elif hs > 0 and 0 <= Azw < 90:
        Ax = calc_Axm(WSSize, Azw, hs)
    else:
        Ax = 0
        
    return Ax


# ## C. 天空日射・反射日射の効果係数 (仕様書6.4, 6.5)

# ### C.1 形態係数算定のための関数$\;f_A\;$ (仕様書6.4 式(23))
# 
# $$ \begin{align}
# f_{A}(x_a,x_b,y_a,y_b,z_a) &= \frac{x_b \sqrt{y_b^2+z_a^2}}{2} \tan^{-1} \frac{x_b}{\sqrt{y_b^2+z_a^2}}
# - \frac{x_b \sqrt{y_a^2+z_a^2}}{2} \tan^{-1} \frac{x_b}{\sqrt{y_a^2+z_a^2}} \\
# &- \frac{x_a \sqrt{y_b^2+z_a^2}}{2} \tan^{-1} \frac{x_a}{\sqrt{y_b^2+z_a^2}}
# + \frac{x_a \sqrt{y_a^2+z_a^2}}{2} \tan^{-1} \frac{x_a}{\sqrt{y_a^2+z_a^2}} \\
# &+ \frac{x_b^2 - y_b^2 - z_a^2}{8} \log (x_b^2 + y_b^2 + z_a^2) - \frac{x_b^2 - y_a^2 - z_a^2}{8} \log (x_b^2 + y_a^2 + z_a^2) \\
# &- \frac{x_a^2 - y_b^2 - z_a^2}{8} \log (x_a^2 + y_b^2 + z_a^2) + \frac{x_a^2 - y_a^2 - z_a^2}{8} \log (x_a^2 + y_a^2 + z_a^2) 
# \qquad \qquad \qquad (23) \\
# \end{align} $$

# In[18]:


""" 式(23) """
import numpy as np
    
def calc_fa_atan(x, y, z):

    if y**2 + z**2 > 0:
        fa_atan = x * ( y**2 + z**2 ) **0.5 / 2 * np.arctan( x / ( y**2 + z**2 ) **0.5 )
    else:
        fa_atan = 0
    
    return fa_atan


def calc_fa_log(x, y, z):

    if x**2 + y**2 + z**2 > 0:
        fa_log = ( x**2 - y**2 - z**2 ) / 8 * np.log( x**2 + y**2 + z**2 )
    else:
        fa_log = 0
        
    return fa_log


def calc_fa(xa, xb, ya, yb, za):
    
    fa = ( calc_fa_atan(xb, yb, za) - calc_fa_atan(xb, ya, za)
         - calc_fa_atan(xa, yb, za) + calc_fa_atan(xa, ya, za)
         + calc_fa_log(xb, yb, za)  - calc_fa_log(xb, ya, za)
         - calc_fa_log(xa, yb, za)  + calc_fa_log(xa, ya, za) )
  
    return fa    


# ### C.2 天空に対する形態係数 (仕様書6.4 式(22))
# 
# $$ \begin{align}
# \phi_{j,y+} = (\pi A_{wind,j})^{-1} \{ & \hspace{2pt} f_A(x_{3y+}, x_2 + x_{3y+}, y_1, y_1 + y_2, z_{y+}) + f_A(y_{1x+}, y_{1x+} + y_2, x_{3}, x_2 + x_{3}, z_{x+}) \\
# + & \hspace{2pt} f_A(x_{1y+}, x_{1y+} + x_2, y_1, y_1 + y_2, z_{y+}) + f_A(y_{1x-}, y_{1x-} + y_2, x_{1}, x_1 + x_2, z_{x-}) \\ 
# + & \hspace{2pt} f_A(x_{3}, x_2 + x_{3}, y_1, y_1 + y_2, 0) + f_A(y_{1}, y_{1} + y_2, x_{3}, x_2 + x_{3}, 0) \\
# + & \hspace{2pt} f_A(x_{1}, x_{1} + x_2, y_1, y_1 + y_2, 0) + f_A(y_{1}, y_{1} + y_2, x_{1}, x_1 + x_2, 0) \\ 
# - & \hspace{2pt} f_A(x_{3y+}, x_2 + x_{3y+}, y_1, y_1 + y_2, 0) - f_A(y_{1x+}, y_{1x+} + y_2, x_{3}, x_2 + x_{3}, 0) \\
# - & \hspace{2pt} f_A(x_{1y+}, x_{1y+} + x_2, y_1, y_1 + y_2, 0) - f_A(y_{1x-}, y_{1x-} + y_2, x_{1}, x_1 + x_2, 0) \qquad (22) \\
# \end{align} $$

# In[19]:


""" 式(22) """
import numpy as np
    
def calc_phiyp(WSSize):
    
    [X1, X2, X3, X1yp, X1ym, X3yp, X3ym, Y1, Y2, Y3, Y1xp, Y1xm, Y3xp, Y3xm, Zxp, Zxm, Zyp, Zym] = WSSize
    
    phiyp = ( 1 / ( np.pi * X2 * Y2 )
            * ( calc_fa(X3yp, X2 + X3yp, Y1, Y1 + Y2, Zyp) + calc_fa(Y1xp, Y1xp + Y2, X3, X2 + X3, Zxp) 
              + calc_fa(X1yp, X1yp + X2, Y1, Y1 + Y2, Zyp) + calc_fa(Y1xm, Y1xm + Y2, X1, X1 + X2, Zxm)  
              + calc_fa(X3,   X2  +  X3, Y1, Y1 + Y2, 0  ) + calc_fa(Y1,   Y1  +  Y2, X3, X2 + X3, 0  ) 
              + calc_fa(X1,   X1  +  X2, Y1, Y1 + Y2, 0  ) + calc_fa(Y1,   Y1  +  Y2, X1, X1 + X2, 0  )     
              - calc_fa(X3yp, X2 + X3yp, Y1, Y1 + Y2, 0  ) - calc_fa(Y1xp, Y1xp + Y2, X3, X2 + X3, 0  ) 
              - calc_fa(X1yp, X1yp + X2, Y1, Y1 + Y2, 0  ) - calc_fa(Y1xm, Y1xm + Y2, X1, X1 + X2, 0  ) ) ) 
    phiyp = max(0, min(phiyp, 0.5))    #負値は0に、0.5を超える場合は0.5で頭打ち
    
    return phiyp     


# ### C.3 天空日射の効果係数 (仕様書6.4 式(21))
# 
# $$\gamma_{isr,j,y+} = 2 \phi_{j,y+} \qquad (21) $$

# In[20]:


""" 式(21) """

def calc_gammayp(WSSize):
    
    gammayp = 2 * calc_phiyp(WSSize)
         
    return gammayp     


# ### C.4 地面に対する形態係数 (仕様書6.5 式(25))
# 
# $$ \begin{align}
# \phi_{j,y-} = (\pi A_{wind,j})^{-1} \{ & \hspace{2pt} f_A(x_{3y-}, x_2 + x_{3y-}, y_3, y_2 + y_3, z_{y-}) + f_A(y_{3x+}, y_2 + y_{3x+}, x_{3}, x_2 + x_{3}, z_{x+}) \\
# + & \hspace{2pt} f_A(x_{1y-}, x_{1y-} + x_2, y_3, y_2 + y_3, z_{y-}) + f_A(y_{3x-}, y_2 + y_{3x-}, x_{1}, x_1 + x_2, z_{x-}) \\ 
# + & \hspace{2pt} f_A(x_{3}, x_2 + x_{3}, y_3, y_2 + y_3, 0) + f_A(y_{3}, y_{2} + y_3, x_{3}, x_2 + x_{3}, 0) \\
# + & \hspace{2pt} f_A(x_{1}, x_{1} + x_2, y_3, y_2 + y_3, 0) + f_A(y_{3}, y_{2} + y_3, x_{1}, x_1 + x_2, 0) \\ 
# - & \hspace{2pt} f_A(x_{3y-}, x_2 + x_{3y-}, y_3, y_2 + y_3, 0) - f_A(y_{3x+}, y_2 + y_{3x+}, x_{3}, x_2 + x_{3}, 0) \\
# - & \hspace{2pt} f_A(x_{1y-}, x_{1y-} + x_2, y_3, y_2 + y_3, 0) - f_A(y_{3x-}, y_2 + y_{3x-}, x_{1}, x_1 + x_2, 0) \qquad (25) \\
# \end{align} $$

# In[21]:


""" 式(25) """
import numpy as np

def calc_phiym(WSSize):
    
    [X1, X2, X3, X1yp, X1ym, X3yp, X3ym, Y1, Y2, Y3, Y1xp, Y1xm, Y3xp, Y3xm, Zxp, Zxm, Zyp, Zym] = WSSize
    
    phiym = ( 1 / ( np.pi * X2 * Y2 )
            * ( calc_fa(X3ym, X2 + X3ym, Y3, Y2 + Y3, Zym) + calc_fa(Y3xp, Y2 + Y3xp, X3, X2 + X3, Zxp) 
              + calc_fa(X1ym, X1ym + X2, Y3, Y2 + Y3, Zym) + calc_fa(Y3xm, Y2 + Y3xm, X1, X1 + X2, Zxm)  
              + calc_fa(X3,   X2  +  X3, Y3, Y2 + Y3, 0  ) + calc_fa(Y3,   Y2  +  Y3, X3, X2 + X3, 0  ) 
              + calc_fa(X1,   X1  +  X2, Y3, Y2 + Y3, 0  ) + calc_fa(Y3,   Y2  +  Y3, X1, X1 + X2, 0  )     
              - calc_fa(X3ym, X2 + X3ym, Y3, Y2 + Y3, 0  ) - calc_fa(Y3xp, Y2 + Y3xp, X3, X2 + X3, 0  ) 
              - calc_fa(X1ym, X1ym + X2, Y3, Y2 + Y3, 0  ) - calc_fa(Y3xm, Y2 + Y3xm, X1, X1 + X2, 0  ) ) ) 
    phiym = max(0, min(phiym, 0.5))    #負値は0に、0.5を超える場合は0.5で頭打ち
    
    return phiym    


# ### C.5 反射日射の効果係数 (仕様書6.5 式(24))
# 
# $$\gamma_{isr,j,y-} = 2 \phi_{j,y-} \qquad (24) $$

# In[22]:


""" 式(24) """

def calc_gammaym(WSSize):

    gammaym = 2 * calc_phiym(WSSize)
         
    return gammaym


# ## D. 地点と日射量
# 
# ### D.1 地点データ読み込み
# 
# - `\Zone.csv`
#   - $1$行目はヘッダ：地域区分, 都市, 緯度, 経度, 日射量ファイル名, 暖房開始日, 暖房終了日, 冷房開始日, 冷房終了日
# 
#   - $2～9$行目：$1～8$地域の「地域区分, 都市, 緯度, 経度, 日射量ファイル名, 暖房開始日, 暖房終了日, 冷房開始日, 冷房終了日」
#     - 暖冷房開始日終了日の書式：$5$桁もしくは$6$桁の数字
#       - 気象データの日時表記と同じ
#       - 後ろから$2$桁：時刻
#       - 後ろから$3～4$桁：日
#       - 後ろから$5～6$桁：月
#       - 気象データの「暖房$1$冷房$2$」の設定は上書き処理される → `\Zone.csv` の設定が優先
# 
# 
# 

# In[23]:


""" 地点データを \Zone.csv から読み込む """
# \地域区分+日射量データ窓面入射角特性.xlsx "地域区分"シート
#   → \SCFConfig01 下の地点データファイル \Zone.csv を作成 → 読み込み
import pandas as pd
import sys

def input_Point(ClimateZone, Path00, FileName00):
    # Path00 = "./SCFConfig01/"
    # FileName00 = "Zone.csv"

    csv_input = pd.read_csv(filepath_or_buffer=Path00+FileName00, encoding="ms932", sep=",")
    if csv_input.columns[0]!="地域区分":
        sys.exit("地点データではありません")
    [Zone, City, Latitude, Longitude, SRFileName, HStart, HEnd, CStart, CEnd]         = ["None", "None", 0, 0, "None", 0, 0, 0, 0]
    for i in range(len(csv_input)):
        if ClimateZone == csv_input.values[i][0] or ClimateZone == csv_input.values[i][1]:
            [Zone, City, Latitude, Longitude, SRFileName, HStart, HEnd, CStart, CEnd] = csv_input.values[i]            

    if Zone=="None":
        sys.exit("地点データがありません")
        
    return [Zone, City, Latitude, Longitude, SRFileName, HStart, HEnd, CStart, CEnd]
        # ここで返るSRFileNameはファイル名のみ


# ### D.2 気象データ読み込み
# 
# - `\SRforSCF_\*\*.csv`
#   - \*\*は、$1～8$地域に対応して$01$～$08$が入る(基本)。
#     - `\Zone.csv` に追加して対応することは可能。   
#     
#     
#   - 法線面直達日射量、水平面天空日射量、暖房期or冷房期の判別タグ(暖房期:$1$, 冷房期:$2$, 非空調期:$0$)
#     - 日射量の単位は$[kcal/(m2h)]$ → 効果係数算定には問題ないのでそのまま使用している
#     
#     
#   - $1$行目はヘッダ：`\SRforSCF_\*\*.csv`(ファイル名), 法線面直達日射量, 水平面天空日射量, 暖房$1$_冷房$2$
#   - $2～8762$行目：日時, 法線面直達日射量, 水平面天空日射量, 暖房$1$_冷房$2$
#     - $2$行目が$1$月$1$日$0$時、$8762$行目が$12$月$31$日$24$時。$1$時間間隔。全$8761$データ
#     - 気象データファイル中の「暖房$1$冷房$2$」の設定は、`\Zone.csv` の設定で上書きされる → 現時点で意味なし
#     - $1$列目の日時から「月」を算出
# 

# In[24]:


""" 気象データを \SRforSCF_**.csv から読み込む """
# \地域区分+日射量データ窓面入射角特性.xlsx "SRforSCF_**.csv"シート 
#   → \SCFConfig01 下の地点データファイル \SRforSCF_**.csv を作成 → 読み込み
import pandas as pd
import sys

def input_SRData(Path00, FileName00, HStart, HEnd, CStart, CEnd):
    # Path00 = "./SCFConfig01/"
    # FileName00 = "SRforSCF_**.csv"

    csv_input = pd.read_csv(filepath_or_buffer=Path00+FileName00, encoding="ms932", sep=",")
    if csv_input.columns[0]!=FileName00[-len(csv_input.columns[0]):]:
        sys.exit("データが違います")

    """ \Zone.csv の設定で、暖房期,冷房期,非空調期を割り当て """
    # 元ファイルの4列目はなかったことになる。
    [HeatingPeriod, CoolingPeriod, NonACPeriod] = [1, 2, 0]
    for i in range(len(csv_input)):
        if CStart <= csv_input.values[i][0] <= CEnd:
            csv_input.values[i][3] = CoolingPeriod
        elif CStart > CEnd and CStart <= csv_input.values[i][0]:
            csv_input.values[i][3] = CoolingPeriod
        elif CStart > CEnd and csv_input.values[i][0] <= CEnd:
            csv_input.values[i][3] = CoolingPeriod
        elif HStart <= csv_input.values[i][0] <= HEnd:
            csv_input.values[i][3] = HeatingPeriod        
        elif HStart > HEnd and HStart <= csv_input.values[i][0]:
            csv_input.values[i][3] = HeatingPeriod
        elif HStart > HEnd and csv_input.values[i][0] <= HEnd:
            csv_input.values[i][3] = HeatingPeriod
        else:
            csv_input.values[i][3] = NonACPeriod        
    
    return csv_input
        # csv_input は SRHour に渡される

    
""" 「月」の計算 """
def calc_Month(MMDDTT):
    
    Month = MMDDTT // 10000
    
    return Month


# ### D.3 正時±30分で太陽が地平線上にある時間刻み数のカウント (仕様書6.2 式(3)及び図5中の$n_H$の計算)
# 
# - 算定ツール標準の時間分割数$n_{\Delta t}$は、$6$ ($10$分刻み)
#   - $1$時間を$1$分割もしくは$2$以上の偶数で分割する
#   
#   
# - 正時$\pm 30$分間で太陽が地平線上にある(太陽高度$>0$)時間刻み数をカウントして、$n_H$を計算

# In[25]:


""" 式(3),図5中のNhの計算 """
import sys

def calc_Nh(Latitude, Longitude, NDay, NHour, NDT):

    deltad = calc_deltad(NDay)
    eed = calc_eed(NDay)
   
    Nh = 0
    if NDT == 1:
        if calc_sinh(Latitude, deltad, calc_Tdt(Longitude, eed, NHour)) > 0:
            Nh += 1        
    elif NDT > 0 and NDT % 2 == 0:
        # 1/1 0時, 各日24時において、赤緯と均時差にズレが生じるが、
        # 白夜でなければ日が昇らないので実質影響なし → 放置
        sinh0 = [calc_sinh(Latitude, deltad, calc_Tdt(Longitude, eed, NHour + m / NDT)) 
                                                  for m in range(-int(NDT/2),int(NDT/2)+1)]
        Nh = ( sum(x > 0 for x in sinh0) 
             - (0.5 if sinh0[0] > 0 else 0) - (0.5 if sinh0[int(NDT)] > 0 else 0) )
        
        #以下のコメントアウト部分を上の式二つで置き換え    
        #if calc_sinh(Latitude, deltad, calc_Tdt(Longitude, eed, NHour - 0.5)) > 0: #正時の30分前        
        #    Nh += 0.5
        #if calc_sinh(Latitude, deltad, calc_Tdt(Longitude, eed, NHour + 0.5)) > 0: #正時の30分後
        #    Nh += 0.5
        #for m in range(int(NDT/2+1), int(NDT)):
        #    if calc_sinh(Latitude, deltad, calc_Tdt(Longitude, eed, NHour -1 + m / NDT)) > 0:
        #        Nh += 1
        #for m in range(0, int(NDT/2)):
        #    if calc_sinh(Latitude, deltad, calc_Tdt(Longitude, eed, NHour + m / NDT)) > 0:
        #        Nh += 1            
    else:
        sys.exit("1時間あたりの時間分割数は1もしくは2以上の偶数とする必要があります")

    return Nh


# ### D.4 $1/n_{\Delta t}$時間間隔での日射量 (仕様書6.2 式(3)の計算, 図5参照)
# 
# - 算定ツール標準の時間分割数$n_{\Delta t}$は、$6$ ($10$分刻み)
# - 正時$\pm 30$分間で太陽が地平線上にある(太陽高度$>0$)時間刻み数をカウントして、$n_H$を計算
# - 法線面直達日射量, 水平面天空日射量をそれぞれに適用し、時間刻みにおける日射量を算定する

# In[26]:


""" 式(3)の S'HM の計算 """
def calc_Sdhm(MM, NDT, sinh, Sh, Shp, Nh, Nhp):
    
    Sdhm = 0
    if sinh > 0:
        Sdhm = ( ( (Sh  / Nh  if MM <= NDT/2 and Nh  > 0 else 0)
                 + (Shp / Nhp if MM >= NDT/2 and Nhp > 0 else 0) )
               / (2 if MM == NDT/2 else 1) )
        #以下のコメントアウト部分を上式で置き換え
        #if MM < NDT/2 and Nh > 0:
        #    Sdhm = Sh / Nh
        #elif MM == NDT/2:
        #    if Nh > 0:
        #        Sdhm += Sh / Nh / 2
        #    if Nhp > 0:
        #        Sdhm += Shp / Nhp / 2            
        #elif NDT/2 < MM and Nhp > 0:
        #    Sdhm = Shp / Nhp
    else:
        Sdhm = 0
        
    return Sdhm  


# ### D.5 窓ガラスの入射角特性読み込み
# 
# - `\IncidentAngleCharacteristics.csv`
#   - $1$行目はヘッダ："入射角特性", $\eta_{max}$, $\eta_{isr}$, $k_0$～$k_7$
#   - $2$行目以下に入射角特性のデータを記入
#     - $1$列：$ID$
#     - $2$列：直達日射に対する入射角特性最大値(入射角$0$) $\eta_{max}$
#     - $3$列：天空・反射日射に対する入射角特性(遮蔽なしの場合) $\eta_{isr}$
#     - $4～11$列：$\eta_{j,d,t}(\theta_{j,d,t})$ 算定式の係数$k_n$$(n=0～7)$
#        $$\eta_{j,d,t}(\theta_{j,d,t}) = \sum_{n=0}^7 k_n \cos^n \theta_{j,d,t}$$
#   - デフォルトとして、以下を設定
#     - $ID=0$：日よけ効果係数内で入射角特性非考慮 → $\eta_{j,d,t}(\theta_{j,d,t}) = 1$
#     - $ID=1$：解説書の入射角特性(「平成25年度省エネルギー基準に準拠した算定・判断の方法及び解説 I 非住宅建築物 第二版(連合印刷センター, 平成26年○月○日)」, pp.168-170, 式(2.1.25),(2.1.28),(2.1.32))
#        $$\eta_{j,d,t}(\theta_{j,d,t}) = 2.3920 \cos \theta_{j,d,t} -3.8636 \cos^3 \theta_{j,d,t} + 3.7568 \cos^5 \theta_{j,d,t} - 1.3952 \cos^7 \theta_{j,d,t} $$
#     - 他の特性を入れる際には、$ID$を違えて、`\IncidentAngleCharacteristics.csv` に追加する。

# In[27]:


""" 入射角特性データセットを \IncidentAngleCharacteristics.csv から読み込む  """
# \地域区分+日射量データ+窓面入射角特性.xlsx "入射角特性"シート
#   → \IncidentAngleCharacteristics.csv を作成 → 読み込み
import pandas as pd
import sys

def input_IncidentAngleCharacteristics(ID, Path00, FileName00):
    # Path00 = "./SCFConfig01/"
    # FileName00 = "IncidentAngleCharacteristics.csv"

    csv_input = pd.read_csv(filepath_or_buffer=Path00+FileName00, encoding="ms932", sep=",")
    if csv_input.columns[0]!="入射角特性":
        sys.exit("ファイル内に貼り付けたテスト条件が違います")

    ID0="none"    
    for i in range(len(csv_input)):
        if ID == csv_input.values[i][0]:
            [ID0, etamax, etaisr, etakk] = [csv_input.values[i][0],csv_input.values[i][1]
                                            ,csv_input.values[i][2],csv_input.values[i][3:11]]
    if ID0=="none":
        sys.exit("指定したIDの入射角特性がありません")            
    
    return [ID0, etamax, etaisr, etakk]


# ### D.6 直達日射に対する窓ガラスの入射角特性 (緑本非住宅第二版pp.169 式(2.1.28)準拠)
# 
# - 「平成25年度省エネルギー基準に準拠した算定・判断の方法及び解説 I 非住宅建築物 第二版(連合印刷センター, 平成26年○月○日)」pp.168-170参照
# - 入射角$\theta_{j,d,t}[deg]$, 太陽高度$h_{S,d,t}[deg]$, 太陽方位角$A_{ZS,d,t}[deg]$, 外壁$j$の方位角$A_{ZW,i}[deg]$
# - 入射角特性算定式の係数$k_n$$(n=0～7)$, 日付$d$時刻$t$における入射角特性$\eta_{j,d,t}$
# 
# $$\cos \theta_{j,d,t} = \cos h_{S,d,t} \cos (A_{ZS,d,t} - A_{ZW,i}) \qquad (2.1.26) $$
# $$\eta_{j,d,t}(\theta_{j,d,t}) = \sum_{n=0}^7 k_n \cos^n \theta_{j,d,t} \qquad (2.1.28') $$

# In[28]:


""" 式(2.1.26),(2.1.28') """
import numpy as np

def calc_costheta(Azwjdt, cosh):

    costheta = max(cosh * np.cos(np.radians(Azwjdt)),0)

    return costheta


def calc_etajdt(costheta, etakk):

    etajdt = sum([etakk[i]*costheta**i for i in range(len(etakk))])
        
    return etajdt


# ## E. 日よけ効果係数
# 
# ### E.1 日よけ効果係数の算出 (仕様書6.1 式(2)準拠)
# 
# - $\gamma_{dsr,j,x+}$：窓等$j$に対して太陽が$x+$側に位置する際の直達日射に対する日よけ効果係数(式$(13)$)
# - $\gamma_{dsr,j,x-}$：窓等$j$に対して太陽が$x-$側に位置する際の直達日射に対する日よけ効果係数(式$(17)$)
# - $\gamma_{isr,j,y+}$：天空日射に対する日よけ効果係数(式$(21)$) → C.3参照
# - $\gamma_{isr,j,y-}$：反射日射に対する日よけ効果係数(式$(24)$) → C.5参照
# - $\eta_{j,d,t}$：日付$d$、時刻$t$における窓等$j$の入射角特性 → D.6参照
# - $S_{D,d,t}$：日付$d$、時刻$t$における法線面直達日射量$[kcal/m^2]$
# - $S_{S,d,t}$：日付$d$、時刻$t$における水平面天空日射量$[kcal/m^2]$
# - $\eta_{isr}$：天空・反射日射に対する入射角特性(遮蔽なしの場合)$[-]$
#   - 緑本、仕様書では$\;\eta_{isr} = 0.808\;$
#   
#   
# - 式(2)のうち、定数$0.5$は垂直面からみた天空・地表面の形態係数、$0.1$は地表面における日射反射率である。
# 
# $$ \begin{align}
# \gamma_{wind,j} 
# = \big[&\hspace{4pt}\gamma_{dsr,j,x+} \times \sum_{-90<A_{ZW,j,d,t}<0, \; h_{s,d,t}>0} 
#         (S_{D,d,t} \; \eta_{j,d,t} \cos h_{s,d,t} \cos A_{ZW,j,d,t}) \\
#   + &\hspace{4pt} \gamma_{dsr,j,x-} \times \sum_{0 \leq A_{ZW,j,d,t}<90, \; h_{s,d,t}>0} 
#         (S_{D,d,t} \; \eta_{j,d,t} \cos h_{s,d,t} \cos A_{ZW,j,d,t}) \\
#   + &\hspace{4pt}\gamma_{isr,j,y+} \times \sum (\eta_{isr} \times 0.5 \times S_{S,d,t} \\
#   + &\hspace{4pt}\gamma_{isr,j,y-} \times \sum \{\eta_{isr} \times 0.1 \times 0.5 
#         \times (S_{S,d,t} + S_{D,d,t} \sin h_{s,d,t}) \} \big] \\
# \big/ \hspace{4pt} \big[&\sum_{-90<A_{ZW,j,d,t}<90, \; h_{s,d,t}>0} (S_{D,d,t} \; \eta_{j,d,t} \cos h_{s,d,t} \cos A_{ZW,j,d,t}) \\
#     + &\hspace{4pt}\sum \{\eta_{isr} \times 0.5 \times S_{S,d,t} 
#         +  \eta_{isr} \times 0.1 \times 0.5 \times (S_{S,d,t} + S_{D,d,t} \sin h_{s,d,t}) \} \big] \qquad (2)
# \\
# \end{align} $$
# 
# - $A_{wind,j}$：窓等$j$の面積($=x_2y_2$)$[m^2]$
# - $A_{wind,j,x+,d,t}$：太陽が$x+$側に位置する日付$d$、時刻$t$において窓等$j$の直達日射が当たる部分の面積$[m^2]$
# - $A_{wind,j,x-,d,t}$：太陽が$x-$側に位置する日付$d$、時刻$t$において窓等$j$の直達日射が当たる部分の面積$[m^2]$
# 
# $$ \begin{align}
# \gamma_{dsr,j,x+} 
# = \big\{ & \sum_{-90<A_{ZW,j,d,t}<0, \; h_{s,d,t}>0} (A_{wind,j,x+,d,t}\hspace{2pt}  
#         S_{D,d,t} \; \eta_{j,d,t} \cos h_{s,d,t} \cos A_{ZW,j,d,t}) \big\} \\
# \big/ \big\{ &\hspace{2pt} A_{wind,j} \sum_{-90<A_{ZW,j,d,t}<0, \; h_{s,d,t}>0} (S_{D,d,t} \; \eta_{j,d,t} \cos h_{s,d,t} \cos A_{ZW,j,d,t})  \} \big\} \qquad (13)
# \\
# \end{align} $$
# 
# $$ \begin{align}
# \gamma_{dsr,j,x-} 
# = \big\{ & \sum_{0 \leq A_{ZW,j,d,t}<90, \; h_{s,d,t}>0} (A_{wind,j,x-,d,t} \hspace{2pt} 
#         S_{D,d,t} \; \eta_{j,d,t} \cos h_{s,d,t} \cos A_{ZW,j,d,t}) \big\} \\
# \big/ \big\{ &\hspace{2pt} A_{wind,j} \sum_{0 \leq A_{ZW,j,d,t}<90, \; h_{s,d,t}>0} (S_{D,d,t} \; \eta_{j,d,t} \cos h_{s,d,t} \cos A_{ZW,j,d,t})  \} \big\} \qquad (17)
# \\
# \end{align} $$
# 

# 
# ### E.2 式(2)分子分母中の時刻$t$における各成分の計算

# In[29]:


""" 直達成分の分子分母 """
def calc_dendsr00(Sddhm, etajdt, costheta, Awj): #(2)式分母中直達成分×窓面積
    
    dendsr00 = Awj * Sddhm * etajdt * costheta

    return dendsr00


def calc_numdsr00(Sddhm, etajdt, costheta, Ax): #(2)式分子中直達成分×窓面積←(13)(17)式

    numdsr00 = Ax * Sddhm * etajdt * costheta
    
    return numdsr00


""" 天空成分の分子分母 """
def calc_denisryp00(Ssdhm, etaisr, Awj): #(2)式分母中天空成分×窓面積
    
    denisryp00 = Awj * etaisr * 0.5 * Ssdhm
    
    return denisryp00


def calc_numisryp00(Ssdhm, etaisr, Awj, gammayp): #(2)式分子中天空成分×窓面積

    numisryp00 = gammayp * Awj * etaisr * 0.5 * Ssdhm
    # numisryp00 = gammayp * denisryp00
    
    return numisryp00    


""" 反射成分の分子分母 """
def calc_denisrym00(Sddhm, Ssdhm, etaisr, sinh, Awj): #(2)式分母中反射成分×窓面積
    
    denisrym00 = Awj * etaisr * 0.1 * 0.5 * (Ssdhm + Sddhm * sinh)
    
    return denisrym00


def calc_numisrym00(Sddhm, Ssdhm, etaisr, sinh, Awj, gammaym): #(2)式分子中反射成分×窓面積

    numisrym00 = gammaym * Awj * etaisr * 0.1 * 0.5 * (Ssdhm + Sddhm * sinh)
    # numisrym00 = gammaym * denisrym00
    
    return numisrym00


""" 分子分母それぞれの和 """
def calc_SCF00(Sddhm, Ssdhm, etajdt, etaisr, costheta, sinh, Awj, Ax, gammayp, gammaym):
    
    SCF00 = [0,0]
    
    # SCF00[0]:(2)式の分子への加算分(窓面積をかけた値として)
    SCF00[0] = ( calc_numdsr00(Sddhm, etajdt, costheta, Ax)
               + calc_numisryp00(Ssdhm, etaisr, Awj, gammayp)
               + calc_numisrym00(Sddhm, Ssdhm, etaisr, sinh, Awj, gammaym) )
    # SCF00[1]:(2)式の分母への加算分(窓面積をかけた値として)
    SCF00[1] = ( calc_dendsr00(Sddhm, etajdt, costheta, Awj)             
               + calc_denisryp00(Ssdhm, etaisr, Awj)
               + calc_denisrym00(Sddhm, Ssdhm, etaisr, sinh, Awj) )

    return SCF00


# 
# ### E.3 式(2)分子分母の期間積算処理
# 
# - 式(2)中の時刻$t$における分子および分母を期間積算する。
#   - 暖房期、冷房期、非空調期
#   - 各月毎
#   - 各月、各時刻毎
#     - 時刻毎は正時$\pm30$分としたいが、本ツール内でのデータの取り扱いからずれる。前時刻の`MM` $=$ `NDT`$/2～$同時刻の`MM` $=$ `NDT`$/2-1$の間を正時の積算値として計上する。
#     
# 
# - `SCF[h][i][j]`に格納して積算
#   - `h`$=0$：分子, `h`$=1$：分母, `h`$=2$：効果係数
#   - `i`$=-2$：冷房期積算, `i`$=-1$：暖房期積算, `i`$=0$：非空調期積算, `i`$=1～12$：各月で積算
#   - `j`$=-1$：日積算, `j`$=0～24$：各時刻で積算  
# 
#   
# 
# 

# In[30]:


""" 期間積算処理 → 日よけ効果係数算出 """

def Output_ShadingCorrectionFactor(SRHour, NDT, SCF01):

    SCF = [[[0 for j in range(-1,25)] for i in range(-2,13)] for h in range(0,3)] 
        # SCF[h][i][j]に格納して積算
        # h = 0：分子, h = 1：分母, h = 2：効果係数
        # i = -2：冷房期積算, i = -1：暖房期積算, i = 0：非空調期積算,
        #         i = 1～12：各月で積算
        # j = -1：日積算, j = 0～24：各時刻で積算  
        
    """ +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 1時間のループ ++++ """       
    for Hour00 in range(8760):      # ← 12/31 24時は回さないのでHour00=8759がループの最後
        
        """ 日時の計算(A.1) """
        [NDay,NHour] = calc_NDayNHour(Hour00)           
        
        """ 「月」の計算(D.2) """
        Month = calc_Month(SRHour.values[Hour00][0])        
 
        """ +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ (1/NDT)時間のループ ++++ """ 
            # MM:1時間の内の(1/NDT)間隔の順番, 正時がMM=0, MM=0～NDT-1
        for MM in range(NDT):
            """ 時刻TTの計算(A.1) """
                # TT:時間分割MM毎の時刻[hour]
            TT = calc_TT(NHour, NDT, MM)
            
            """ 時刻Hour01の計算(A.1) """
                # Hour01:前時刻のMM=NDT/2～同時刻のMM=NDT/2-1の時間分割を、その正時に属するものとして
                #        扱うための時刻 → 効果係数算定のための積算用
            Hour01 = calc_Hour01(TT)
            
            """ 暖冷房期間のタグ """            
            HCTag = -SRHour.values[Hour00][3]

            """ 分子分母の期間,月,時間毎の積算(E.3) """
            for h in range(0, 2):
                SCF[h][HCTag][Hour01] += SCF01[h][Hour00][MM]
                SCF[h][HCTag][-1]     += SCF01[h][Hour00][MM]        
                SCF[h][Month][Hour01] += SCF01[h][Hour00][MM]
                SCF[h][Month][-1]     += SCF01[h][Hour00][MM]                      
                    # SCF[h][i][j]に格納して積算
                    # h = 0：分子, h = 1：分母, h = 2：効果係数
                    # i = -2：冷房期積算, i = -1：暖房期積算, i = 0：非空調期積算,
                    #         i = 1～12：各月で積算
                    # j = -1：日積算, j = 0～24：各時刻で積算             
                    
    """ 期間,月,時間毎の日よけ効果係数算出(E.1) """
    for i in range(-2,13):
        for j in range(-1,25):
            if SCF[1][i][j] != 0:
                SCF[2][i][j] = SCF[0][i][j] / SCF[1][i][j]
            else:
                SCF[2][i][j] = 0

    return SCF


