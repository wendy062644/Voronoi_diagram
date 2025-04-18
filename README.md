## NSYSU_Design and Analysis of Algorithm (演算法設計與分析 CSE510)

### 介紹
以分治法(Divide and Conquer)為基礎的 Voronoi Diagram 建構演算法，並搭配 Python GUI 介面呈現互動式的視覺化結果。  
核心使用 C++ 計算，視覺化使用 Python (Tkinter) 製作。  

### 架構
C++ (voronoi.cpp)	建構 Voronoi Diagram，處理幾何計算。  
Python (index.py)	Tkinter GUI：提供畫布操作、點擊互動、結果顯示。  
透過 points.txt、lines.txt 作為中介資料交換。  

### 使用方法
pip install numpy scipy  
g++ voronoi.cpp -o voronoi  
python index.py  

### 操作方式
1. 點擊畫布新增點  
2. 按下 Execute Current Points 生成 Voronoi 圖  
3. 使用 Prev/Next Step 查看分治細節  
4. 可選擇匯出結果 (Save Output)  

### 分治法演算法說明
基底情況 (Base case)：
2點：直接建構垂直平分線
3點：計算三邊垂直平分線並求交點。

遞迴步驟 (Recursive Step)：  
將點依據 x 座標排序後，左右分開處理。  
對左右子問題分別遞迴計算 Voronoi Diagram。  

合併步驟 (Merge Step)：  
找到左、右邊的 Convex Hull  
計算 Hyperplane，逐步構建左右子圖的合併邊界。  

### 預期改進方向
引入更多異常情況處理。  
改善 C++ 程式的時間與空間效率。  

[書面報告連結](https://wendy062644.github.io/voronoi_report/index.html "link")
