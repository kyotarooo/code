import pandas as pd

# fileの読み込み及び表示
file = "c:\\m\\bigdata.csv"
df = pd.read_csv(file, index_col=0, encoding='Shift-JIS')
print('==============================')
print('all data(問1):')
print(df)
print('==============================')
print(f'データ行(問2)  :{len(df)}')
print('==============================')

#欠損の処理
df_Nan = df.isnull()
Nan_num = df_Nan.sum()
print(f'欠損値の個数(問3):\n  {Nan_num}')
print('==============================')
print('欠損値のないデータ(問4):\n')
df_dropna = df.dropna()
print(df_dropna)
print('==============================')
Nan_dropna_num = len(df_dropna)
print(f'データ行(問5): {Nan_dropna_num}')
print('==============================')

# data = {'国語':[50,50,None,40], '数学':[80,None,None,50]}
# idx = ['Alice', 'Bob', 'Chie', 'Dan']
# dfA = pd.DataFrame(data,index=idx)
# print(dfA)
# print('\n')
# print(dfA.isnull())
# print('\n')
# print(dfA.isnull().sum())
# dfB = dfA.dropna()
# print(dfB)
# print(dfA.mean())
# dfC = dfA.fillna(dfA.mean())
# print(dfC)

# import pandas as pd
# # ディクショナリによるカラム名も含めたdf作成
# data = {     
# 	'国語' : [90,50,None,40],
# 	'数学' : [80,None,None,50]
# }
# idx =  ['Alice', 'Bob ', 'Chie ', 'Dan']
# dfA = pd.DataFrame(data, index=idx)
# print(dfA)
# print('\n')
# print(dfA.isnull())  #NanならTrueを表示
# print('\n')
# print ( dfA.isnull().sum() )  # sumでTrueの個数を合計


# import pandas as pd

# file ='C:\\m\\' + 'test.csv'    #文字列の連結．絶対パス形式での指定
# print(file)
# dfA = pd.read_csv(file, encoding='UTF-8',index_col=0)   #変数の場合は'  'はナシ
# print(dfA)

# dfB = dfA[dfA['国語']>80]
# print(dfB)
# print('\n')

# dfC = dfA[(dfA['国語']>80)&(dfA['数学']>=80) ]
# print(dfC)


