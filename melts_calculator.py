#!/usr/bin/env python3
"""
MELTS Calculator - Windows Compatible Version
A thermodynamic equilibrium calculator for magmatic systems
Based on MELTS Excel implementation
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
import pandas as pd
import json
import os
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.style as mplstyle
mplstyle.use('seaborn-v0_8-darkgrid')

class MELTSCalculator:
    def __init__(self, root):
        self.root = root
        self.root.title("MELTS Calculator - Windows Edition")
        self.root.geometry("1200x800")
        
        # 스타일 설정
        self.setup_styles()
        
        # 데이터 초기화
        self.init_data()
        
        # UI 생성
        self.create_ui()
        
        # 초기값 설정
        self.set_default_values()
        
    def setup_styles(self):
        """UI 스타일 설정"""
        style = ttk.Style()
        style.theme_use('clam')
        
        # 색상 테마
        self.colors = {
            'bg': '#f0f0f0',
            'fg': '#333333',
            'button': '#4CAF50',
            'button_hover': '#45a049',
            'error': '#f44336',
            'warning': '#ff9800',
            'info': '#2196F3'
        }
        
        # 버튼 스타일
        style.configure('Calculate.TButton', 
                       background=self.colors['button'],
                       foreground='white',
                       borderwidth=0,
                       focuscolor='none',
                       font=('Arial', 10, 'bold'))
        
        style.map('Calculate.TButton',
                 background=[('active', self.colors['button_hover'])])
    
    def init_data(self):
        """데이터 구조 초기화"""
        # 산화물 성분
        self.oxides = [
            'SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 
            'FeO', 'MnO', 'MgO', 'NiO', 'CoO', 
            'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 
            'CO2', 'SO3', 'Cl2O-1', 'F2O-1'
        ]
        
        # 상(phase) 목록
        self.phases = [
            'liquid', 'feldspar', 'water', 'quartz', 'tridymite',
            'cristobalite', 'clinopyroxene', 'orthopyroxene', 
            'olivine', 'spinel', 'rhm-oxide', 'garnet', 'biotite',
            'cordierite', 'hornblende', 'cummingtonite', 'pigeonite'
        ]
        
        # 입력 데이터 저장
        self.composition_vars = {}
        self.phase_vars = {}
        self.results_data = {}
        
    def create_ui(self):
        """메인 UI 생성"""
        # 메인 컨테이너
        main_container = ttk.Frame(self.root, padding="10")
        main_container.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # 탭 생성
        self.notebook = ttk.Notebook(main_container)
        self.notebook.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # 각 탭 생성
        self.create_input_tab()
        self.create_conditions_tab()
        self.create_phases_tab()
        self.create_results_tab()
        self.create_plots_tab()
        
        # 하단 버튼
        button_frame = ttk.Frame(main_container)
        button_frame.grid(row=1, column=0, pady=10, sticky=(tk.W, tk.E))
        
        ttk.Button(button_frame, text="Calculate", 
                  command=self.calculate,
                  style='Calculate.TButton').pack(side=tk.LEFT, padx=5)
        
        ttk.Button(button_frame, text="Reset", 
                  command=self.reset_values).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(button_frame, text="Save Results", 
                  command=self.save_results).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(button_frame, text="Load Data", 
                  command=self.load_data).pack(side=tk.LEFT, padx=5)
        
        # 상태바
        self.status_var = tk.StringVar(value="Ready")
        status_bar = ttk.Label(main_container, textvariable=self.status_var, 
                             relief=tk.SUNKEN, anchor=tk.W)
        status_bar.grid(row=2, column=0, sticky=(tk.W, tk.E))
        
        # 그리드 가중치 설정
        main_container.columnconfigure(0, weight=1)
        main_container.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
    
    def create_input_tab(self):
        """조성 입력 탭"""
        input_frame = ttk.Frame(self.notebook)
        self.notebook.add(input_frame, text="Composition")
        
        # 스크롤 가능한 프레임
        canvas = tk.Canvas(input_frame, bg='white')
        scrollbar = ttk.Scrollbar(input_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # 헤더
        ttk.Label(scrollable_frame, text="Oxide", font=('Arial', 10, 'bold')).grid(
            row=0, column=0, padx=5, pady=5, sticky=tk.W)
        ttk.Label(scrollable_frame, text="Weight %", font=('Arial', 10, 'bold')).grid(
            row=0, column=1, padx=5, pady=5)
        ttk.Label(scrollable_frame, text="Moles", font=('Arial', 10, 'bold')).grid(
            row=0, column=2, padx=5, pady=5)
        
        # 산화물 입력 필드
        for i, oxide in enumerate(self.oxides):
            row = i + 1
            
            # 산화물 이름
            ttk.Label(scrollable_frame, text=oxide).grid(
                row=row, column=0, padx=5, pady=2, sticky=tk.W)
            
            # Weight % 입력
            wt_var = tk.DoubleVar(value=0.0)
            self.composition_vars[oxide] = {
                'wt': wt_var,
                'moles': tk.StringVar(value="0.000")
            }
            
            wt_entry = ttk.Entry(scrollable_frame, textvariable=wt_var, width=10)
            wt_entry.grid(row=row, column=1, padx=5, pady=2)
            wt_entry.bind('<KeyRelease>', lambda e, ox=oxide: self.update_moles(ox))
            
            # Moles 표시 (읽기 전용)
            moles_label = ttk.Label(scrollable_frame, 
                                   textvariable=self.composition_vars[oxide]['moles'])
            moles_label.grid(row=row, column=2, padx=5, pady=2)
        
        # 합계 표시
        ttk.Separator(scrollable_frame, orient='horizontal').grid(
            row=len(self.oxides)+1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        
        ttk.Label(scrollable_frame, text="Total:", font=('Arial', 10, 'bold')).grid(
            row=len(self.oxides)+2, column=0, padx=5, pady=5, sticky=tk.W)
        
        self.total_wt_var = tk.StringVar(value="0.00")
        ttk.Label(scrollable_frame, textvariable=self.total_wt_var, 
                 font=('Arial', 10, 'bold')).grid(
            row=len(self.oxides)+2, column=1, padx=5, pady=5)
        
        # 정규화 버튼
        ttk.Button(scrollable_frame, text="Normalize to 100%", 
                  command=self.normalize_composition).grid(
            row=len(self.oxides)+3, column=0, columnspan=3, pady=10)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
    
    def create_conditions_tab(self):
        """계산 조건 탭"""
        conditions_frame = ttk.Frame(self.notebook)
        self.notebook.add(conditions_frame, text="Conditions")
        
        # 온도-압력 조건
        tp_frame = ttk.LabelFrame(conditions_frame, text="Temperature-Pressure Conditions", 
                                 padding="10")
        tp_frame.grid(row=0, column=0, padx=10, pady=10, sticky=(tk.W, tk.E))
        
        # 온도
        ttk.Label(tp_frame, text="Temperature (°C):").grid(row=0, column=0, sticky=tk.W)
        self.temp_var = tk.DoubleVar(value=800.0)
        temp_entry = ttk.Entry(tp_frame, textvariable=self.temp_var, width=15)
        temp_entry.grid(row=0, column=1, padx=5)
        
        # 압력
        ttk.Label(tp_frame, text="Pressure (bar):").grid(row=1, column=0, sticky=tk.W)
        self.pressure_var = tk.DoubleVar(value=200.0)
        pressure_entry = ttk.Entry(tp_frame, textvariable=self.pressure_var, width=15)
        pressure_entry.grid(row=1, column=1, padx=5)
        
        # 산소 분압
        ttk.Label(tp_frame, text="log fO₂:").grid(row=2, column=0, sticky=tk.W)
        self.logfo2_var = tk.DoubleVar(value=0.0)
        logfo2_entry = ttk.Entry(tp_frame, textvariable=self.logfo2_var, width=15)
        logfo2_entry.grid(row=2, column=1, padx=5)
        
        # 버퍼 선택
        ttk.Label(tp_frame, text="fO₂ Buffer:").grid(row=2, column=2, sticky=tk.W, padx=(20,0))
        self.buffer_var = tk.StringVar(value="None")
        buffer_combo = ttk.Combobox(tp_frame, textvariable=self.buffer_var, 
                                   values=["None", "QFM", "NNO", "IW", "MW", "FMQ"],
                                   width=10)
        buffer_combo.grid(row=2, column=3, padx=5)
        
        # 계산 모드
        mode_frame = ttk.LabelFrame(conditions_frame, text="Calculation Mode", padding="10")
        mode_frame.grid(row=1, column=0, padx=10, pady=10, sticky=(tk.W, tk.E))
        
        self.calc_mode = tk.StringVar(value="equilibrium")
        ttk.Radiobutton(mode_frame, text="Equilibrium", 
                       variable=self.calc_mode, 
                       value="equilibrium").grid(row=0, column=0, padx=5)
        ttk.Radiobutton(mode_frame, text="Fractional Crystallization", 
                       variable=self.calc_mode, 
                       value="fractional").grid(row=0, column=1, padx=5)
        
        # 시퀀스 설정
        seq_frame = ttk.LabelFrame(conditions_frame, text="T-P Sequence", padding="10")
        seq_frame.grid(row=2, column=0, padx=10, pady=10, sticky=(tk.W, tk.E))
        
        ttk.Label(seq_frame, text="Initial T (°C):").grid(row=0, column=0, sticky=tk.W)
        self.t_start_var = tk.DoubleVar(value=1200.0)
        ttk.Entry(seq_frame, textvariable=self.t_start_var, width=10).grid(row=0, column=1)
        
        ttk.Label(seq_frame, text="Final T (°C):").grid(row=0, column=2, sticky=tk.W, padx=(20,0))
        self.t_end_var = tk.DoubleVar(value=700.0)
        ttk.Entry(seq_frame, textvariable=self.t_end_var, width=10).grid(row=0, column=3)
        
        ttk.Label(seq_frame, text="Step (°C):").grid(row=0, column=4, sticky=tk.W, padx=(20,0))
        self.t_step_var = tk.DoubleVar(value=10.0)
        ttk.Entry(seq_frame, textvariable=self.t_step_var, width=10).grid(row=0, column=5)
        
        self.sequence_enable = tk.BooleanVar(value=False)
        ttk.Checkbutton(seq_frame, text="Enable Sequence", 
                       variable=self.sequence_enable).grid(row=1, column=0, columnspan=6, pady=5)
    
    def create_phases_tab(self):
        """상(Phase) 선택 탭"""
        phases_frame = ttk.Frame(self.notebook)
        self.notebook.add(phases_frame, text="Phases")
        
        info_label = ttk.Label(phases_frame, 
                             text="Select phases to include in calculations:",
                             font=('Arial', 10, 'italic'))
        info_label.grid(row=0, column=0, columnspan=4, pady=10)
        
        # 상 체크박스들
        for i, phase in enumerate(self.phases):
            row = (i // 4) + 1
            col = i % 4
            
            var = tk.BooleanVar(value=True if phase in ['liquid', 'feldspar', 'quartz'] else False)
            self.phase_vars[phase] = var
            
            cb = ttk.Checkbutton(phases_frame, text=phase, variable=var)
            cb.grid(row=row, column=col, padx=10, pady=5, sticky=tk.W)
        
        # 전체 선택/해제 버튼
        button_frame = ttk.Frame(phases_frame)
        button_frame.grid(row=(len(self.phases)//4)+2, column=0, columnspan=4, pady=20)
        
        ttk.Button(button_frame, text="Select All", 
                  command=lambda: self.toggle_all_phases(True)).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Deselect All", 
                  command=lambda: self.toggle_all_phases(False)).pack(side=tk.LEFT, padx=5)
    
    def create_results_tab(self):
        """결과 표시 탭"""
        results_frame = ttk.Frame(self.notebook)
        self.notebook.add(results_frame, text="Results")
        
        # 결과 트리뷰
        columns = ('Property', 'System', 'Liquid', 'Solids', 'Water')
        self.results_tree = ttk.Treeview(results_frame, columns=columns, show='tree headings')
        
        # 컬럼 설정
        self.results_tree.heading('#0', text='Phase')
        for col in columns:
            self.results_tree.heading(col, text=col)
            self.results_tree.column(col, width=120)
        
        # 스크롤바
        v_scroll = ttk.Scrollbar(results_frame, orient='vertical', 
                               command=self.results_tree.yview)
        h_scroll = ttk.Scrollbar(results_frame, orient='horizontal', 
                               command=self.results_tree.xview)
        self.results_tree.configure(yscrollcommand=v_scroll.set, 
                                  xscrollcommand=h_scroll.set)
        
        # 배치
        self.results_tree.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        v_scroll.grid(row=0, column=1, sticky=(tk.N, tk.S))
        h_scroll.grid(row=1, column=0, sticky=(tk.W, tk.E))
        
        results_frame.columnconfigure(0, weight=1)
        results_frame.rowconfigure(0, weight=1)
        
        # 결과 요약 프레임
        summary_frame = ttk.LabelFrame(results_frame, text="Summary", padding="10")
        summary_frame.grid(row=2, column=0, columnspan=2, pady=10, sticky=(tk.W, tk.E))
        
        self.summary_text = tk.Text(summary_frame, height=6, width=80)
        self.summary_text.pack(fill=tk.BOTH, expand=True)
    
    def create_plots_tab(self):
        """그래프 탭"""
        plots_frame = ttk.Frame(self.notebook)
        self.notebook.add(plots_frame, text="Plots")
        
        # 그래프 선택
        control_frame = ttk.Frame(plots_frame)
        control_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=10)
        
        ttk.Label(control_frame, text="Plot Type:").pack(side=tk.LEFT, padx=5)
        self.plot_type = tk.StringVar(value="phase_diagram")
        plot_combo = ttk.Combobox(control_frame, textvariable=self.plot_type,
                                 values=["phase_diagram", "composition", "properties"],
                                 width=20)
        plot_combo.pack(side=tk.LEFT, padx=5)
        
        ttk.Button(control_frame, text="Update Plot", 
                  command=self.update_plot).pack(side=tk.LEFT, padx=10)
        
        # 그래프 캔버스
        self.fig, self.ax = plt.subplots(figsize=(8, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=plots_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def update_moles(self, oxide):
        """몰수 업데이트"""
        try:
            wt = self.composition_vars[oxide]['wt'].get()
            # 분자량 (간단한 예시 값)
            mol_weights = {
                'SiO2': 60.08, 'TiO2': 79.87, 'Al2O3': 101.96,
                'Fe2O3': 159.69, 'Cr2O3': 151.99, 'FeO': 71.84,
                'MnO': 70.94, 'MgO': 40.30, 'NiO': 74.69,
                'CoO': 74.93, 'CaO': 56.08, 'Na2O': 61.98,
                'K2O': 94.20, 'P2O5': 141.94, 'H2O': 18.02,
                'CO2': 44.01, 'SO3': 80.06, 'Cl2O-1': 70.90,
                'F2O-1': 37.00
            }
            
            if oxide in mol_weights:
                moles = wt / mol_weights[oxide]
                self.composition_vars[oxide]['moles'].set(f"{moles:.3f}")
            
            # 총합 업데이트
            total = sum(self.composition_vars[ox]['wt'].get() for ox in self.oxides)
            self.total_wt_var.set(f"{total:.2f}")
            
        except:
            pass
    
    def normalize_composition(self):
        """조성을 100%로 정규화"""
        try:
            total = sum(self.composition_vars[ox]['wt'].get() for ox in self.oxides)
            if total > 0:
                for oxide in self.oxides:
                    current = self.composition_vars[oxide]['wt'].get()
                    normalized = (current / total) * 100
                    self.composition_vars[oxide]['wt'].set(normalized)
                    self.update_moles(oxide)
                
                self.status_var.set("Composition normalized to 100%")
        except Exception as e:
            messagebox.showerror("Error", f"Normalization failed: {str(e)}")
    
    def toggle_all_phases(self, select):
        """모든 상 선택/해제"""
        for phase_var in self.phase_vars.values():
            phase_var.set(select)
    
    def calculate(self):
        """MELTS 계산 수행"""
        try:
            self.status_var.set("Calculating...")
            self.root.update()
            
            # 입력 데이터 수집
            composition = {ox: self.composition_vars[ox]['wt'].get() 
                          for ox in self.oxides}
            
            temperature = self.temp_var.get()
            pressure = self.pressure_var.get()
            logfo2 = self.logfo2_var.get()
            
            selected_phases = [phase for phase, var in self.phase_vars.items() if var.get()]
            
            # 시퀀스 계산 여부 확인
            if self.sequence_enable.get():
                self.calculate_sequence(composition, selected_phases)
            else:
                # 단일 점 계산
                results = self.perform_calculation(composition, temperature, 
                                                 pressure, logfo2, selected_phases)
                self.display_results(results)
            
            self.status_var.set("Calculation completed")
            
        except Exception as e:
            messagebox.showerror("Calculation Error", str(e))
            self.status_var.set("Calculation failed")
    
    def perform_calculation(self, composition, T, P, logfo2, phases):
        """실제 MELTS 계산 수행 (간단한 시뮬레이션)"""
        # 이것은 실제 MELTS 알고리즘의 간단한 근사입니다
        # 실제 구현에서는 복잡한 열역학 계산이 필요합니다
        
        results = {
            'temperature': T,
            'pressure': P,
            'logfo2': logfo2,
            'system': {},
            'phases': {}
        }
        
        # 시스템 전체 속성 계산
        total_mass = sum(composition.values())
        results['system']['mass'] = total_mass
        results['system']['density'] = 2.65  # 근사값
        results['system']['G'] = -1700 - T * 0.1  # 간단한 근사
        results['system']['H'] = -1400 - T * 0.08
        results['system']['S'] = 270 + T * 0.01
        results['system']['Cp'] = 140 + T * 0.005
        results['system']['V'] = total_mass / results['system']['density']
        
        # 액체 상 계산 (간단한 예시)
        if 'liquid' in phases:
            liquid_frac = min(1.0, max(0.0, (T - 700) / 500))  # 온도에 따른 액체 분율
            results['phases']['liquid'] = {
                'mass': total_mass * liquid_frac,
                'composition': composition.copy(),
                'density': 2.4,
                'viscosity': 10 ** (4 - T/300)  # 간단한 점성도 모델
            }
        
        # 고체 상 계산
        solid_mass = total_mass * (1 - liquid_frac) if 'liquid' in phases else total_mass
        
        # 주요 광물 상 분배 (간단한 예시)
        if solid_mass > 0:
            if 'feldspar' in phases and composition.get('Al2O3', 0) > 10:
                results['phases']['feldspar'] = {
                    'mass': solid_mass * 0.3,
                    'formula': 'NaAlSi3O8',
                    'density': 2.62
                }
            
            if 'quartz' in phases and composition.get('SiO2', 0) > 60:
                results['phases']['quartz'] = {
                    'mass': solid_mass * 0.2,
                    'formula': 'SiO2',
                    'density': 2.65
                }
            
            if 'orthopyroxene' in phases and composition.get('MgO', 0) > 0:
                results['phases']['orthopyroxene'] = {
                    'mass': solid_mass * 0.15,
                    'formula': 'MgSiO3',
                    'density': 3.2
                }
        
        return results
    
    def calculate_sequence(self, composition, phases):
        """온도 시퀀스 계산"""
        try:
            t_start = self.t_start_var.get()
            t_end = self.t_end_var.get()
            t_step = self.t_step_var.get()
            pressure = self.pressure_var.get()
            logfo2 = self.logfo2_var.get()
            
            temperatures = np.arange(t_start, t_end - t_step/2, -abs(t_step))
            sequence_results = []
            
            # 진행 상황 창
            progress_window = tk.Toplevel(self.root)
            progress_window.title("Calculating Sequence...")
            progress_window.geometry("300x100")
            
            progress_var = tk.DoubleVar()
            progress_bar = ttk.Progressbar(progress_window, variable=progress_var, 
                                         maximum=len(temperatures))
            progress_bar.pack(padx=20, pady=20, fill=tk.X)
            
            status_label = ttk.Label(progress_window, text="Processing...")
            status_label.pack()
            
            for i, T in enumerate(temperatures):
                results = self.perform_calculation(composition, T, pressure, 
                                                 logfo2, phases)
                sequence_results.append(results)
                
                progress_var.set(i + 1)
                status_label.config(text=f"Temperature: {T:.1f}°C")
                progress_window.update()
            
            progress_window.destroy()
            
            # 시퀀스 결과 플롯
            self.plot_sequence_results(sequence_results)
            
        except Exception as e:
            messagebox.showerror("Sequence Calculation Error", str(e))
    
    def display_results(self, results):
        """결과를 트리뷰에 표시"""
        # 기존 결과 삭제
        for item in self.results_tree.get_children():
            self.results_tree.delete(item)
        
        # 시스템 속성
        system_item = self.results_tree.insert('', 'end', text='System Properties')
        props = ['mass', 'density', 'G', 'H', 'S', 'Cp', 'V']
        for prop in props:
            if prop in results['system']:
                values = ['', prop, f"{results['system'][prop]:.3f}", '', '']
                self.results_tree.insert(system_item, 'end', values=values)
        
        # 상별 결과
        phases_item = self.results_tree.insert('', 'end', text='Phase Properties')
        for phase, data in results['phases'].items():
            phase_item = self.results_tree.insert(phases_item, 'end', text=phase)
            for key, value in data.items():
                if isinstance(value, (int, float)):
                    values = ['', key, '', f"{value:.3f}", '']
                    self.results_tree.insert(phase_item, 'end', values=values)
        
        # 요약 텍스트 업데이트
        self.update_summary(results)
        
        # 트리 확장
        self.results_tree.item(system_item, open=True)
        self.results_tree.item(phases_item, open=True)
    
    def update_summary(self, results):
        """결과 요약 업데이트"""
        self.summary_text.delete(1.0, tk.END)
        
        summary = f"Temperature: {results['temperature']}°C\n"
        summary += f"Pressure: {results['pressure']} bar\n"
        summary += f"log fO₂: {results['logfo2']}\n\n"
        
        if 'liquid' in results['phases']:
            liquid_mass = results['phases']['liquid']['mass']
            total_mass = results['system']['mass']
            liquid_frac = (liquid_mass / total_mass) * 100
            summary += f"Liquid fraction: {liquid_frac:.1f}%\n"
            summary += f"Solid fraction: {100-liquid_frac:.1f}%\n"
        
        self.summary_text.insert(1.0, summary)
    
    def update_plot(self):
        """그래프 업데이트"""
        plot_type = self.plot_type.get()
        
        self.ax.clear()
        
        if plot_type == "phase_diagram":
            self.plot_phase_diagram()
        elif plot_type == "composition":
            self.plot_composition()
        elif plot_type == "properties":
            self.plot_properties()
        
        self.canvas.draw()
    
    def plot_phase_diagram(self):
        """상 다이어그램 그리기"""
        # 예시 데이터
        T = np.linspace(700, 1200, 100)
        liquid = 100 * (T - 700) / 500
        solid = 100 - liquid
        
        self.ax.fill_between(T, 0, liquid, alpha=0.5, label='Liquid')
        self.ax.fill_between(T, liquid, 100, alpha=0.5, label='Solid')
        
        self.ax.set_xlabel('Temperature (°C)')
        self.ax.set_ylabel('Phase Fraction (%)')
        self.ax.set_title('Phase Diagram')
        self.ax.legend()
        self.ax.grid(True, alpha=0.3)
    
    def plot_composition(self):
        """조성 막대 그래프"""
        oxides = []
        values = []
        
        for oxide in self.oxides:
            wt = self.composition_vars[oxide]['wt'].get()
            if wt > 0:
                oxides.append(oxide)
                values.append(wt)
        
        if oxides:
            self.ax.bar(oxides, values)
            self.ax.set_xlabel('Oxide')
            self.ax.set_ylabel('Weight %')
            self.ax.set_title('Composition')
            self.ax.tick_params(axis='x', rotation=45)
    
    def plot_properties(self):
        """물성 그래프"""
        # 예시: 온도에 따른 점성도
        T = np.linspace(700, 1200, 100)
        log_visc = 4 - T/300
        
        self.ax.plot(T, log_visc)
        self.ax.set_xlabel('Temperature (°C)')
        self.ax.set_ylabel('log Viscosity (Pa·s)')
        self.ax.set_title('Melt Viscosity vs Temperature')
        self.ax.grid(True, alpha=0.3)
    
    def plot_sequence_results(self, results):
        """시퀀스 결과 플롯"""
        self.notebook.select(4)  # Plots 탭으로 이동
        
        self.ax.clear()
        
        temperatures = [r['temperature'] for r in results]
        liquid_fracs = []
        
        for r in results:
            if 'liquid' in r['phases']:
                lf = r['phases']['liquid']['mass'] / r['system']['mass'] * 100
            else:
                lf = 0
            liquid_fracs.append(lf)
        
        self.ax.plot(temperatures, liquid_fracs, 'b-', linewidth=2)
        self.ax.set_xlabel('Temperature (°C)')
        self.ax.set_ylabel('Liquid Fraction (%)')
        self.ax.set_title('Crystallization Path')
        self.ax.grid(True, alpha=0.3)
        self.ax.invert_xaxis()
        
        self.canvas.draw()
    
    def set_default_values(self):
        """기본값 설정"""
        # 기본 조성 (화강암질)
        default_comp = {
            'SiO2': 74.39,
            'TiO2': 0.18,
            'Al2O3': 13.55,
            'Fe2O3': 0.36,
            'FeO': 0.98,
            'MgO': 0.50,
            'CaO': 1.43,
            'Na2O': 3.36,
            'K2O': 5.09,
            'H2O': 5.04
        }
        
        for oxide, value in default_comp.items():
            if oxide in self.composition_vars:
                self.composition_vars[oxide]['wt'].set(value)
                self.update_moles(oxide)
    
    def reset_values(self):
        """모든 값 초기화"""
        response = messagebox.askyesno("Reset", "Reset all values to defaults?")
        if response:
            # 조성 초기화
            for oxide in self.oxides:
                self.composition_vars[oxide]['wt'].set(0.0)
                self.update_moles(oxide)
            
            # 조건 초기화
            self.temp_var.set(800.0)
            self.pressure_var.set(200.0)
            self.logfo2_var.set(0.0)
            
            # 기본값 설정
            self.set_default_values()
            
            self.status_var.set("Values reset to defaults")
    
    def save_results(self):
        """결과 저장"""
        filename = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                data = {
                    'timestamp': datetime.now().isoformat(),
                    'composition': {ox: self.composition_vars[ox]['wt'].get() 
                                  for ox in self.oxides},
                    'conditions': {
                        'temperature': self.temp_var.get(),
                        'pressure': self.pressure_var.get(),
                        'logfo2': self.logfo2_var.get()
                    },
                    'phases': {phase: var.get() 
                             for phase, var in self.phase_vars.items()},
                    'results': self.results_data
                }
                
                with open(filename, 'w') as f:
                    json.dump(data, f, indent=2)
                
                self.status_var.set(f"Results saved to {filename}")
                messagebox.showinfo("Success", "Results saved successfully!")
                
            except Exception as e:
                messagebox.showerror("Save Error", str(e))
    
    def load_data(self):
        """데이터 불러오기"""
        filename = filedialog.askopenfilename(
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                with open(filename, 'r') as f:
                    data = json.load(f)
                
                # 조성 불러오기
                if 'composition' in data:
                    for oxide, value in data['composition'].items():
                        if oxide in self.composition_vars:
                            self.composition_vars[oxide]['wt'].set(value)
                            self.update_moles(oxide)
                
                # 조건 불러오기
                if 'conditions' in data:
                    self.temp_var.set(data['conditions'].get('temperature', 800))
                    self.pressure_var.set(data['conditions'].get('pressure', 200))
                    self.logfo2_var.set(data['conditions'].get('logfo2', 0))
                
                # 상 설정 불러오기
                if 'phases' in data:
                    for phase, value in data['phases'].items():
                        if phase in self.phase_vars:
                            self.phase_vars[phase].set(value)
                
                self.status_var.set(f"Data loaded from {filename}")
                messagebox.showinfo("Success", "Data loaded successfully!")
                
            except Exception as e:
                messagebox.showerror("Load Error", str(e))


def main():
    """메인 함수"""
    root = tk.Tk()
    app = MELTSCalculator(root)
    root.mainloop()


if __name__ == "__main__":
    main()