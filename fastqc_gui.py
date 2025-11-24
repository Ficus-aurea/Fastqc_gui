import tkinter as tk
from tkinter import filedialog, messagebox
import ttkbootstrap as ttk
from ttkbootstrap.constants import *
from ttkbootstrap.widgets.scrolled import ScrolledText
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import os
from pathlib import Path


BG_COLOR = "#121212"       
FG_COLOR = "#FFD1A9"        
ACCENT_COLOR = "#FF8C42"   
PLOT_BG = "#1E1E1E"         
PLOT_FG = "#FFD1A9"         
FONT_MAIN = ("Calibri", 11)
FONT_HEADER = ("Calibri", 12, "bold")
FONT_MONO = ("Consolas", 10) 

class FastqReader:
    def __init__(self, path): self.path = path
    def __enter__(self): return self
    def __exit__(self, exc_type, exc_val, exc_tb): pass
    def read(self): 
        import random
        for _ in range(100):
            seq = "".join(random.choices("ATGC", k=random.randint(100, 150)))
            qual = [random.randint(20, 40) for _ in range(len(seq))]
            yield type('obj', (object,), {'sequence': seq, 'quality': qual})

class FastQCApp:

    def __init__(self, root):
        self.root = root
        self.root.title("FastQC Pro (Dark Mode)")
        self.root.geometry("1000x800")
        
        
        self.style = ttk.Style()
        self.style.configure('.', font=FONT_MAIN, background=BG_COLOR, foreground=FG_COLOR)
        self.style.configure('TFrame', background=BG_COLOR)
        self.style.configure('TLabel', background=BG_COLOR, foreground=FG_COLOR, font=FONT_MAIN)
        self.style.configure('TButton', font=FONT_HEADER)
        self.style.configure('TNotebook', background=BG_COLOR)
        self.style.configure('TNotebook.Tab', font=FONT_MAIN)
        
        main_container = ttk.Frame(self.root, padding=20)
        main_container.pack(fill=BOTH, expand=YES)

        top_frame = ttk.Frame(main_container)
        top_frame.pack(side=TOP, fill=X, pady=(0, 20))

        self.btn_open = ttk.Button(
            top_frame, 
            text=" ・❥・ Открыть FASTQ файл", 
            command=self.load_file, 
            bootstyle="warning", 
            width=25
        )
        self.btn_open.pack(side=LEFT, padx=(0, 15))

        self.lbl_file = ttk.Label(
            top_frame, 
            text="Файл не выбран", 
            font=("Calibri", 11, "italic"),
            foreground=FG_COLOR
        )
        self.lbl_file.pack(side=LEFT, pady=5)

        self.notebook = ttk.Notebook(main_container, bootstyle="dark")
        self.notebook.pack(fill=BOTH, expand=YES)

       
        self.tab_summary = ttk.Frame(self.notebook, padding=10)
        self.tab_len = ttk.Frame(self.notebook, padding=10)
        self.tab_qual = ttk.Frame(self.notebook, padding=10)
        self.tab_content = ttk.Frame(self.notebook, padding=10)

        self.notebook.add(self.tab_summary, text=" ✧.* Сводка ")
        self.notebook.add(self.tab_len, text=" ✧.* Длины ")
        self.notebook.add(self.tab_qual, text=" ✧.* Качество ")
        self.notebook.add(self.tab_content, text=" ✧.* Нуклеотиды ")

        self.text_summary = ScrolledText(
            self.tab_summary, 
            font=("Calibri", 12), 
            padding=10,
            autohide=True,
            bootstyle="dark" 
        )
        
        self.text_summary.text.configure(bg=PLOT_BG, fg=FG_COLOR, insertbackground=FG_COLOR)
        self.text_summary.pack(fill=BOTH, expand=YES)

    def load_file(self):
        file_path = filedialog.askopenfilename(
            filetypes=[("FASTQ files", "*.fastq"), ("Gzip FASTQ", "*.fastq.gz"), ("All files", "*.*")]
        )
        if file_path:
            filename = os.path.basename(file_path)
            self.lbl_file.config(text=f"Файл: {filename}")
            self.run_analysis(file_path)

    def run_analysis(self, file_path):
        try:
            self.text_summary.text.delete(1.0, END) 
            
            for tab in [self.tab_len, self.tab_qual, self.tab_content]:
                for widget in tab.winfo_children():
                    widget.destroy()

            file_path = Path(file_path)
            sequence_lengths = []
            quality_per_position = {}
            base_content_per_position = {}
            total_sequences = 0
            total_length = 0

            self.root.config(cursor="watch")
            self.root.update()

            with FastqReader(file_path) as reader:
                for record in reader.read():
                    seq = record.sequence
                    qual = record.quality
                    seq_len = len(seq)
                    total_sequences += 1
                    total_length += seq_len
                    sequence_lengths.append(seq_len)
                    
                    for i, q in enumerate(qual):
                        if i not in quality_per_position: quality_per_position[i] = []
                        quality_per_position[i].append(q)
                    
                    for i, base in enumerate(seq):
                        if i not in base_content_per_position: base_content_per_position[i] = {"A": 0, "T": 0, "G": 0, "C": 0}
                        if base in "ATGC": base_content_per_position[i][base] += 1

            self.root.config(cursor="")
            
            if total_sequences == 0:
                self.text_summary.text.insert(END, "Файл пуст.")
                return

            mean_length = total_length / total_sequences
            summary_text = (
                f"Анализ завершен для: {file_path.name}\n"
                f"{'-'*40}\n"
                f"Всего последовательностей: {total_sequences}\n"
                f"Всего нуклеотидов:       {total_length}\n"
                f"Средняя длина:           {mean_length:.1f} bp\n"
                f"Минимальная длина:       {min(sequence_lengths)} bp\n"
                f"Максимальная длина:      {max(sequence_lengths)} bp\n"
            )
            self.text_summary.text.insert(END, summary_text)

            self.draw_length_dist(sequence_lengths)
            self.draw_quality_plot(quality_per_position)
            self.draw_content_plot(base_content_per_position)

        except Exception as e:
            self.root.config(cursor="")
            messagebox.showerror("Ошибка", f"Ошибка: {e}")

    def embed_plot(self, figure, parent_frame):
        figure.tight_layout() 
        canvas = FigureCanvasTkAgg(figure, master=parent_frame)
        canvas.draw()
        widget = canvas.get_tk_widget()
        widget.pack(fill=BOTH, expand=YES, padx=5, pady=5)
        return canvas

    def _create_figure(self):
       
        fig = Figure(figsize=(8, 5), dpi=100, facecolor=PLOT_BG) 
        ax = fig.add_subplot(111)
        
        
        ax.set_facecolor(PLOT_BG)
        
       
        ax.spines['bottom'].set_color(FG_COLOR)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color(FG_COLOR)
        
       
        ax.tick_params(axis='x', colors=FG_COLOR)
        ax.tick_params(axis='y', colors=FG_COLOR)
        ax.yaxis.label.set_color(FG_COLOR)
        ax.xaxis.label.set_color(FG_COLOR)
        ax.title.set_color(FG_COLOR)
        
        
        return fig, ax


    def draw_length_dist(self, lengths):
        fig, ax = self._create_figure()
       
        ax.hist(lengths, bins=min(50, len(set(lengths))), color=ACCENT_COLOR, edgecolor=PLOT_BG, alpha=0.8)
        ax.set_title("Распределение длин", fontsize=12, pad=15)
        ax.set_xlabel("Длина (bp)")
        ax.grid(True, linestyle=':', alpha=0.3, color=FG_COLOR)
        self.embed_plot(fig, self.tab_len)

    def draw_quality_plot(self, quality_data):
        fig, ax = self._create_figure()
        if quality_data:
            positions = sorted(quality_data.keys())
            mean_qualities = [sum(quality_data[p]) / len(quality_data[p]) for p in positions]
            
            
            ax.plot(positions, mean_qualities, color=FG_COLOR, linewidth=2)
           
           
            ax.fill_between(positions, 0, 20, color='#e74c3c', alpha=0.2) 
            ax.fill_between(positions, 20, 28, color='#f1c40f', alpha=0.2) 
            ax.fill_between(positions, 28, 42, color='#2ecc71', alpha=0.2) 
            
            ax.set_title("Среднее качество (Phred)", fontsize=12, pad=15)
            ax.set_ylim(0, 42)
        self.embed_plot(fig, self.tab_qual)

    def draw_content_plot(self, content_data):
        fig, ax = self._create_figure()
        if content_data:
            positions = sorted(content_data.keys())
            
            
            colors = {'A': '#46f0f0', 'T': '#ff5f5f', 'G': '#f0f046', 'C': '#5fafff'}
            
            for base, color in colors.items():
                vals = []
                for p in positions:
                    total = sum(content_data[p].values())
                    vals.append((content_data[p][base] / total * 100) if total else 0)
                ax.plot(positions, vals, label=base, color=color, linewidth=1.5)
            
           
            legend = ax.legend(loc='upper right', frameon=False)
            plt.setp(legend.get_texts(), color=FG_COLOR)
            
            ax.set_title("Содержание нуклеотидов %", fontsize=12, pad=15)
            ax.set_ylim(0, 100)
            ax.grid(True, linestyle=':', alpha=0.3, color=FG_COLOR)
        self.embed_plot(fig, self.tab_content)

if __name__ == "__main__":
  
    app = ttk.Window(themename="darkly") 
    FastQCApp(app)
    app.mainloop()
