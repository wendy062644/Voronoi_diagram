import tkinter as tk
from tkinter import filedialog, messagebox
from scipy.spatial import ConvexHull
from tkinter import ttk
import math
import subprocess
import os
import numpy as np

canvas_height = 602  # 定義畫布的高度
file_loaded = False
output_file = False
current_step = 0
hyperplane_files = []
convexhull_files = []
test_cases = []

def clear_file_content(directory):
    try:
        for filename in os.listdir(directory):
            file_path = os.path.join(directory, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)
    except Exception as e:
        messagebox.showerror("Error", f"無法清空資料夾內容: {directory}\n{str(e)}")

def clear_canvas():
    """清空畫布並重置點資料。"""
    canvas.delete("all")
    clear_file_content("hyperplane")
    clear_file_content("convexhull")
    points.clear()

def click_event(event):
    """在畫布上繪製點並記錄座標。"""
    x, y = event.x, event.y
    canvas.create_oval(x-3, y-3, x+3, y+3, fill='red', outline='')
    points.append((x, y))

def draw_voronoi():
    """繪製當前點的Voronoi圖，並標記由外部程序生成的點。"""
    if len(points) < 2:
        messagebox.showwarning("Warning", "需要至少兩個點來繪製Voronoi圖。")
        return
    
    # 將點寫入文件
    with open("points.txt", "w") as file:
        for point in points:
            file.write(f"{point[0]} {canvas_height - point[1]}\n")  # 坐標轉換
    
    # 呼叫C++執行Voronoi並捕捉輸出
    result = subprocess.run(["./voronoi"], capture_output=True, text=True)
    
    load_step_by_step_files()
    
    if result.returncode != 0:
        messagebox.showerror("Error", "Voronoi計算失敗。")
        clear_canvas()
        return

    # 處理 C++ 程序的標準輸出 (假設輸出為 "x y" 格式的點)
    voronoi_output = result.stdout.strip()
    if voronoi_output:
        for line in voronoi_output.splitlines():
            try:
                x, y = map(float, line.split())
                y_canvas = canvas_height - y  # Y軸轉換
                canvas.create_oval(x-3, y_canvas-3, x+3, y_canvas+3, fill='green', outline='')  # 標記輸出點
            except ValueError:
                continue  # 忽略非正常數據行

    # 處理 lines.txt 來畫線
    if os.path.exists("lines.txt"):
        with open("lines.txt", "r") as file:
            for line in file:
                parts = line.split()
                if len(parts) == 4:
                    try:
                        x1, y1, x2, y2 = map(float, parts)
                        y1 = canvas_height - y1  # 坐標轉換
                        y2 = canvas_height - y2
                        canvas.create_line(x1, y1, x2, y2, fill='blue', width=2)
                    except ValueError:
                        continue  # 忽略非正常數據行'

    if os.path.exists("points.txt"):
        with open("points.txt", "r") as file:
            for line in file:
                parts = line.split()
                if len(parts) == 2:
                    try:
                        x, y = map(float, parts)
                        y = canvas_height - y  # 坐標轉換
                        canvas.create_oval(x-3, y-3, x+3, y+3, fill='red', outline='')
                    except ValueError:
                        continue  # 忽略非正常數據行'
    draw_hyperplane()

def load_file():
    """從檔案載入測試資料並存入測資集。"""
    global test_cases, current_case_index, file_loaded  # 確保是全域變數
    file_path = filedialog.askopenfilename(filetypes=[("Text Files", "*.txt")])
    if not file_path:
        return
    
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    test_cases = []
    temp_points = []

    for line in lines:
        line = line.strip()
        if line.startswith("#") or line == "":
            continue

        try:
            num_points = int(line)
            if temp_points:
                test_cases.append(temp_points)
                temp_points = []
            if num_points == 0:
                break
        except ValueError:
            x, y = map(int, line.split())
            temp_points.append((x, y))
    
    if temp_points:
        test_cases.append(temp_points)
    
    current_case_index = 0
    file_loaded = True  # 載入成功
    load_test_case()
    
def load_test_case():
    """載入當前測資並繪製在畫布上。"""
    global current_case_index, test_cases
    clear_canvas()
    if 0 <= current_case_index < len(test_cases):
        for x, y in test_cases[current_case_index]:
            y = canvas_height - y  # 坐標轉換
            canvas.create_oval(x-3, y-3, x+3, y+3, fill='red', outline='')
            points.append((x, y))
        draw_voronoi()
        case_label.config(text=f"Case：{current_case_index + 1} / {len(test_cases)}")
    else:
        messagebox.showinfo("Info", "沒有更多測試資料。")

def next_test_case():
    """顯示下一筆測資。"""
    global current_case_index, current_step
    if output_file:
        if current_case_index < len(draw_test_cases) - 1:
            current_case_index += 1
            load_draw_test_case()
        
    elif current_case_index < len(test_cases) - 1:
        current_case_index += 1
        current_step = 0
        load_test_case()

def previous_test_case():
    """顯示上一筆測資。"""
    global current_case_index, current_step
    if output_file:
        if current_case_index > 0:
            current_case_index -= 1
            load_draw_test_case()
    elif current_case_index > 0:
        current_case_index -= 1
        current_step = 0
        load_test_case()

def execute_current_points():
    """執行當前畫布上的點並生成Voronoi圖。"""
    global current_case_index, test_cases
    canvas.delete("all")
    
    if len(points) < 2:
        messagebox.showwarning("Warning", "需要至少兩個點來繪製Voronoi圖。")
        return
    
    draw_voronoi()
    
    test_cases = []
    test_cases.append(points)
    current_case_index = 0
        
def clear_file_data():
    """清除上傳的測試資料記錄並切換至手動模式。"""
    global test_cases, current_case_index, file_loaded, current_step, hyperplane_files, convexhull_files
    clear_canvas()
    test_cases = []
    hyperplane_files = []
    convexhull_files = []
    current_case_index = -1
    current_step = 0
    file_loaded = False  # 確保手動模式
    case_label.config(text="Case：0 / 0")
    step_label.config(text="Step：0 / 0")

def save_output():
    """將所有測資的點座標和線段輸出至同一檔案。"""
    output_file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text Files", "*.txt")])
    if not output_file_path:
        return
    
    all_results = []  # 存放所有測資結果

    # 如果畫布上有手動繪製的點且未載入檔案，則輸出手繪測資
    if points and not file_loaded:
        all_results.append(f"# 測資")
        process_and_save_current_points(all_results)
    
    # 處理從檔案載入的測資
    if file_loaded:
        for case_index, test_case in enumerate(test_cases):
            all_results.append(f"# 測資 {case_index + 1}")
            process_and_save_case(test_case, all_results)

    # 將所有結果寫入輸出檔案
    with open(output_file_path, "w") as output_file:
        output_file.write("\n".join(all_results))
    
    messagebox.showinfo("Info", "所有測資結果已成功輸出。")

def process_and_save_case(test_case, all_results):
    """處理單筆測資並存入 all_results。"""
    with open("points.txt", "w") as file:
        for x, y in test_case:
            file.write(f"{x} {y}\n")
    
    subprocess.run(["./voronoi"], capture_output=True, text=True)
    
    # 儲存結果
    save_edges_and_points(test_case, all_results)

def process_and_save_current_points(all_results):
    """處理手動測資並存入 all_results。"""
    with open("points.txt", "w") as file:
        for x, y in points:
            file.write(f"{x} {y}\n")
    
    subprocess.run(["./voronoi"], capture_output=True, text=True)
    
    save_edges_and_points(points, all_results)

def save_edges_and_points(points, all_results):
    """儲存邊與點資訊到 all_results。"""
    unique_points = sorted(set(points))
    edges = []
    
    if os.path.exists("lines.txt"):
        with open("lines.txt", "r") as line_file:
            for line in line_file:
                parts = line.split()
                if len(parts) == 4:
                    x1, y1, x2, y2 = map(float, parts)
                    if not any(math.isnan(coord) for coord in (x1, y1, x2, y2)):
                        edges.append((int(x1), int(y1), int(x2), int(y2)))
    
    for point in unique_points:
        all_results.append(f"P {point[0]} {point[1]}")
    for edge in edges:
        all_results.append(f"E {edge[0]} {edge[1]} {edge[2]} {edge[3]}")
    all_results.append("")  # 空行區隔

def draw_from_data(data):
    """解析並根據輸入資料繪製點與線。"""
    for line in data.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):  # 忽略空行和註解
            continue
        
        parts = line.split()
        if parts[0] == "P" and len(parts) == 3:
            try:
                _, x, y = parts
                x, y = int(x), int(y)
                y = canvas_height - y  # Y軸轉換
                canvas.create_oval(x-3, y-3, x+3, y+3, fill='red', outline='')  # 繪製點
            except ValueError:
                continue
        elif parts[0] == "E" and len(parts) == 5:
            try:
                _, x1, y1, x2, y2 = parts
                x1, y1, x2, y2 = map(int, (x1, y1, x2, y2))
                y1, y2 = canvas_height - y1, canvas_height - y2  # Y軸轉換
                canvas.create_line(x1, y1, x2, y2, fill='blue', width=2)  # 繪製線
            except ValueError:
                continue

def load_draw_test_case():
    """載入並繪製 draw_test_cases 的測資。"""
    clear_canvas()
    if 0 <= current_case_index < len(draw_test_cases):
        draw_from_data(draw_test_cases[current_case_index])
        case_label.config(text=f"Case：{current_case_index + 1} / {len(draw_test_cases)}")
    else:
        messagebox.showinfo("Info", "沒有更多測試資料。")

def load_and_draw():
    """從檔案讀取測試資料，將資料分割成多筆測資，並繪製第一筆。"""
    global draw_test_cases, current_case_index, file_loaded, output_file
    file_path = filedialog.askopenfilename(filetypes=[("Text Files", "*.txt")])
    if not file_path:
        return
    try:
        with open(file_path, "r", encoding="utf-8") as file:
            data = file.read()
    except UnicodeDecodeError:
        with open(file_path, "r", encoding="cp950", errors="replace") as file:
            data = file.read()
    
    # 根據 # 分割為多筆測資
    draw_test_cases = [case.strip() for case in data.split("#") if case.strip()]
    current_case_index = 0  # 初始化為第一筆測資
    file_loaded = False  # 避免與檔案測資混淆
    output_file = True  # 載入成功
    load_draw_test_case()

# 讀取 hyperplane 資料夾中的所有檔案
def load_step_by_step_files():
    global hyperplane_files, convexhull_files
    hyperplane_dir = "hyperplane"
    convexhull_dir = "convexhull"
    try:
        hyperplane_files = [f for f in os.listdir(hyperplane_dir) if f.startswith("hyperplane") and f.endswith(".txt")]
        convexhull_files = [f for f in os.listdir(convexhull_dir) if f.startswith("convexhull") and f.endswith(".txt")]
        hyperplane_files.sort()  # 排序檔案以保證順序
        convexhull_files.sort()  # 排序檔案以保證順序
        step_label.config(text=f"Step：{current_step + 1} / {len(hyperplane_files)}")
    except IOError:
        messagebox.showerror("Error", "無法讀取超平面或凸包檔案。")

def draw_hyperplane():
    
    if current_step < 0 or current_step >= len(hyperplane_files):
        return
    
    hyperplane_points = []
    
    step_label.config(text=f"Step：{current_step + 1} / {len(hyperplane_files)}")

    filename = os.path.join("hyperplane", hyperplane_files[current_step])
    try:
        with open(filename, "r") as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) == 2:
                    try:
                        x, y = map(float, parts)
                        hyperplane_points.append((x, canvas_height - y))  # Y軸轉換
                    except ValueError:
                        continue  # 忽略非正常數據行
    except IOError:
        messagebox.showerror("Error", f"無法讀取檔案: {filename}")
        return

    # 用紅色線連接所有hyperplane點
    if len(hyperplane_points) >= 2:
        for i in range(len(hyperplane_points) - 1):
            x1, y1 = hyperplane_points[i]
            x2, y2 = hyperplane_points[i + 1]
            canvas.create_line(x1, y1, x2, y2, fill='red', width=2)
    
    r_points = []
    l_points = []

    filename = os.path.join("convexhull", convexhull_files[current_step])
    try:
        with open(filename, "r") as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) == 3:
                    label, x, y = parts
                    try:
                        x, y = int(x), int(y)
                        if label == "R":
                            r_points.append((x, canvas_height - y))  # Y軸轉換
                        elif label == "L":
                            l_points.append((x, canvas_height - y))
                    except ValueError:
                        continue  # 忽略非正常數據行
    except IOError:
        messagebox.showerror("Error", f"無法讀取檔案: {filename}")
        return

    # 畫出R類型點的凸包
    if len(r_points) >= 3:
        r_points_np = np.array(r_points)
        hull = ConvexHull(r_points_np)
        for simplex in hull.simplices:
            x1, y1 = r_points_np[simplex[0]]
            x2, y2 = r_points_np[simplex[1]]
            canvas.create_line(x1, y1, x2, y2, fill='green', width=2)
    elif len(r_points) == 2:
        x1, y1 = r_points[0]
        x2, y2 = r_points[1]
        canvas.create_line(x1, y1, x2, y2, fill='green', width=2)

    # 畫出L類型點的凸包
    if len(l_points) >= 3:
        l_points_np = np.array(l_points)
        hull = ConvexHull(l_points_np)
        for simplex in hull.simplices:
            x1, y1 = l_points_np[simplex[0]]
            x2, y2 = l_points_np[simplex[1]]
            canvas.create_line(x1, y1, x2, y2, fill='green', width=2)
    elif len(l_points) == 2:
        x1, y1 = l_points[0]
        x2, y2 = l_points[1]
        canvas.create_line(x1, y1, x2, y2, fill='green', width=2)

    # 在畫布上標記所有點
    for x, y in r_points + l_points:
        canvas.create_oval(x - 3, y - 3, x + 3, y + 3, fill='orange', outline='')

def draw_step(step):
    global current_step
    # 確保步驟在有效範圍內
    if step < 0 or step >= len(hyperplane_files):
        return
    current_step = step
    canvas.delete("all")
    draw_voronoi()

def previous_step():
    if current_step > 0:
        draw_step(current_step - 1)

def next_step():
    if current_step < len(hyperplane_files) - 1:
        draw_step(current_step + 1)

# GUI設定
root = tk.Tk()
root.title("Voronoi Diagram")
root.configure(bg='#f0f0f0')

main_frame = ttk.Frame(root, padding="10")
main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

canvas = tk.Canvas(main_frame, width=600, height=canvas_height, bg='white', relief='sunken', bd=2)
canvas.grid(row=0, column=0, columnspan=5, pady=10)

points = []
test_cases = []
current_case_index = -1

# 按鈕配置
clear_button = ttk.Button(main_frame, text="Clear Canvas", command=clear_file_data)
clear_button.grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)

load_button = ttk.Button(main_frame, text="Load Points from File", command=load_file)
load_button.grid(row=1, column=1, padx=5, pady=5, sticky=tk.W)

execute_button = ttk.Button(main_frame, text="Execute Current Points", command=execute_current_points)
execute_button.grid(row=1, column=2, padx=5, pady=5, sticky=tk.W)

load_draw_button = ttk.Button(main_frame, text="Load Output", command=load_and_draw)
load_draw_button.grid(row=1, column=3, padx=5, pady=5, sticky=tk.W)

prev_button = ttk.Button(main_frame, text="Prev Test Case", command=previous_test_case)
prev_button.grid(row=2, column=0, padx=5, pady=5, sticky=tk.W)

next_button = ttk.Button(main_frame, text="Next Test Case", command=next_test_case)
next_button.grid(row=2, column=2, padx=5, pady=5, sticky=tk.W)

save_button = ttk.Button(main_frame, text="Save Output", command=save_output)
save_button.grid(row=2, column=3, padx=5, pady=5, sticky=tk.W)

case_label = ttk.Label(main_frame, text="Case：0 / 0")
case_label.grid(row=2, column=1, padx=5, pady=5, sticky=tk.W)

previous_button = ttk.Button(main_frame, text="Prev Step", command=previous_step)
previous_button.grid(row=3, column=0, padx=10, pady=10)

step_label = ttk.Label(main_frame, text="Step: 0 / 0")
step_label.grid(row=3, column=1, padx=10, pady=10)

next_button = ttk.Button(main_frame, text="Next Step", command=next_step)
next_button.grid(row=3, column=2, padx=10, pady=10)

canvas.bind("<Button-1>", click_event)

root.mainloop()
