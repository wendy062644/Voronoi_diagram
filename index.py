import matplotlib.pyplot as plt

# 假設這是 C++ 輸出的數據
cpp_output = """
E 104.581395 0.000000 0.000000 214.142857
E 290.053603 214.389208 0.000000 458.754177
E 600.000000 405.452055 290.053603 214.389208
E 377.918033 0.000000 290.053603 214.389208
"""

# 解析 C++ 的輸出，提取線段數據
lines = cpp_output.strip().split('\n')
coordinates = []
for line in lines:
    parts = line.split()
    if parts[0] == 'E':  # 只處理以 'E' 開頭的行
        x1, y1, x2, y2 = map(float, parts[1:])
        coordinates.append(((x1, y1), (x2, y2)))

# 打印檢查解析到的座標
print("Parsed Coordinates:")
for (x1, y1), (x2, y2) in coordinates:
    print(f"Start: ({x1}, {y1}), End: ({x2}, {y2})")

# 設置畫布大小
plt.figure(figsize=(8, 8))

# 繪製每個線段
for (x1, y1), (x2, y2) in coordinates:
    # 繪製線段
    plt.plot([x1, x2], [y1, y2], 'b-')  # 'b-' 代表藍色實線
    #plt.scatter(x1, y1, color='red', label='Extra Point')
    #plt.scatter(x2, y2, color='red', label='Extra Point')

# 定義要標記的額外點 (x, y)
#extra_points = [
#    (24, 42), (432, 43), (63, 24), (43, 542), (412, 324),
#    (2, 1), (543, 432), (99, 18), (59, 321), (243, 85),
#    (4, 9), (25, 97), (11, 567), (469, 413), (197, 328)
#]

extra_points = [
    (567, 234), (79, 34), (34, 90), (432, 453), (77, 111)
]

mid_points = [(55.5, 100.5), (255.5, 243.5), (499.5, 343.5), (323, 134)]

# 標記額外的點
for (x, y) in extra_points:
    plt.scatter(x, y, color='orange', label='Extra Point')  # 用橙色標記額外的點

for (x, y) in mid_points:
    plt.scatter(x, y, color='blue', label='Extra Point')  # 用橙色標記額外的點

# 設置坐標範圍，確保視窗大小為 600x600
plt.xlim(0, 600)
plt.ylim(0, 600)

# 添加網格和標題
plt.grid(True)
plt.title('Parsed Lines with Extra Points')

# 顯示圖形
plt.show()
