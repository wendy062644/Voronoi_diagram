# Python Code for GUI (Tkinter for simplicity)
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
import subprocess
import os

def clear_canvas():
    canvas.delete("all")
    points.clear()

def click_event(event):
    x, y = event.x, event.y
    canvas.create_oval(x-3, y-3, x+3, y+3, fill='red', outline='')
    points.append((x, y))

def draw_voronoi():
    if len(points) < 2:
        return
    # Save points to a file
    with open("points.txt", "w") as file:
        for point in points:
            file.write(f"{point[0]} {point[1]}\n")
    
    # Call C++ code for Voronoi computation
    subprocess.run(["./voronoi"])
    
    # Read and draw lines from output file
    if os.path.exists("lines.txt"):
        with open("lines.txt", "r") as file:
            for line in file:
                x1, y1, x2, y2 = map(float, line.split())
                canvas.create_line(x1, y1, x2, y2, fill='blue', width=2)

def load_file():
    file_path = filedialog.askopenfilename(filetypes=[("Text Files", "*.txt")])
    if file_path:
        with open(file_path, "r", encoding="utf-8") as file:
            lines = file.readlines()
            num_points = 0
            points_started = False
            for line in lines:
                line = line.strip()
                # Skip lines starting with "#" or empty lines
                if line.startswith("#") or line == "":
                    continue
                # If the line starts with a digit, it indicates the number of points
                if line[0].isdigit() and not points_started:
                    try:
                        num_points = int(line.split()[0])
                        points_started = True
                    except ValueError:
                        continue
                # If points have started, read exactly the specified number of points
                elif points_started and num_points > 0:
                    parts = line.split()
                    if len(parts) == 2:
                        try:
                            x, y = map(int, parts)
                            canvas.create_oval(x-3, y-3, x+3, y+3, fill='red', outline='')
                            points.append((x, y))
                            num_points -= 1
                        except ValueError:
                            continue
                    # Stop reading if all points are read
                    if num_points == 0:
                        break
                # If the line starts with "P", it indicates an unknown number of points
                elif line.startswith("P"):
                    _, x, y = line.split()
                    x, y = int(x), int(y)
                    canvas.create_oval(x-3, y-3, x+3, y+3, fill='red', outline='')
                    points.append((x, y))
                # If the line starts with "E", it represents a line to be drawn
                elif line.startswith("E"):
                    _, x1, y1, x2, y2 = line.split()
                    x1, y1, x2, y2 = map(int, [x1, y1, x2, y2])
                    canvas.create_line(x1, y1, x2, y2, fill='blue', width=2)

def save_output():
    output_file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text Files", "*.txt")])
    if output_file_path:
        with open(output_file_path, "w") as file:
            # Write points
            for point in points:
                file.write(f"P {point[0]} {point[1]}\n")
            # Write edges
            if os.path.exists("lines.txt"):
                with open("lines.txt", "r") as line_file:
                    for line in line_file:
                        x1, y1, x2, y2 = map(float, line.split())
                        file.write(f"E {int(x1)} {int(y1)} {int(x2)} {int(y2)}\n")

# Create the GUI window
root = tk.Tk()
root.title("Voronoi Diagram")
root.configure(bg='#f0f0f0')

main_frame = ttk.Frame(root, padding="10")
main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

canvas = tk.Canvas(main_frame, width=600, height=600, bg='white', relief='sunken', bd=2)
canvas.grid(row=0, column=0, columnspan=4, pady=10)
points = []

# Add buttons to clear canvas, draw Voronoi diagram, load points from file, and save output
clear_button = ttk.Button(main_frame, text="Clear Canvas", command=clear_canvas)
clear_button.grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
draw_button = ttk.Button(main_frame, text="Draw Voronoi", command=draw_voronoi)
draw_button.grid(row=1, column=1, padx=5, pady=5, sticky=tk.W)
load_button = ttk.Button(main_frame, text="Load Points from File", command=load_file)
load_button.grid(row=1, column=2, padx=5, pady=5, sticky=tk.W)
save_button = ttk.Button(main_frame, text="Save Output", command=save_output)
save_button.grid(row=1, column=3, padx=5, pady=5, sticky=tk.W)

# Bind mouse click to add points on canvas
canvas.bind("<Button-1>", click_event)

root.mainloop()
