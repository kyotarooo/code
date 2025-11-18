#matplotlibをインストールすること

try:
    import matplotlib.pyplot as plt#type: ignore
except ModuleNotFoundError:
    print("matplotlibがインストールされていません。以下のコマンドをターミナルで実行してください。")
    print("pip install matplotlib")
    exit()

import random

sample_list = [0,1,2,3,4,5,6,7,8,9]
random_num = random.sample(sample_list, 3)
answer_num = "".join(map(str, random_num))

match_num_log = []
exist_num_log = []
number_of_trials = 0

print(f'{answer_num}(先生確認用)')

while True:
    
    match_num = 0
    exist_num = 0
    number_of_trials += 1

    while True:
        input_num = input("input three_digit number :")

        if len(input_num)>3:
            print('エラー:三桁の数字を入力してください。')
            continue

        elif len(set(input_num)) != 3:
            print("エラー：数字が重複しています。異なる3つの数字を入力してください。")
        
        else:
            break

    for i in range(3):
        if input_num[i] == answer_num[i]:
            match_num += 1

        elif input_num[i] in answer_num:#配列に特定の値が含まれるかの判定
            exist_num += 1

    print(f'match数 = {match_num}')
    print(f'exist数 = {exist_num}')

    match_num_log.append(match_num)
    exist_num_log.append(exist_num)

    if match_num == 3:
        print(f'おめでとうございます 試行回数：{number_of_trials}')
        break


plt.plot(range(1, len(match_num_log)+1),match_num_log,marker='o', label='match_number')
plt.plot(range(1, len(exist_num_log)+1),exist_num_log,marker='s', label='exist_number')
plt.xlabel('times')
plt.ylabel('match_number/exist_number')
plt.title('line graph')
plt.xticks(range(1, len(match_num_log)+1, 1))  # X軸整数
plt.yticks(range(0, 4, 1))                     # Y軸0〜3までの整数（3桁なので）
plt.legend()
#plt.grid(True)
plt.tight_layout()
plt.show()







