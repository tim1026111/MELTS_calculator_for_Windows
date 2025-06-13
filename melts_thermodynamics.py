"""
MELTS Thermodynamics Module
열역학 계산 및 상평형 모델링을 위한 핵심 모듈
"""

import numpy as np
from scipy.optimize import minimize, root
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import math


@dataclass
class Phase:
    """상(Phase) 데이터 클래스"""
    name: str
    mass: float
    composition: Dict[str, float]
    properties: Dict[str, float]
    formula: str = ""


class MELTSThermodynamics:
    """MELTS 열역학 계산 엔진"""
    
    def __init__(self):
        # 산화물 분자량 (g/mol)
        self.molecular_weights = {
            'SiO2': 60.0843, 'TiO2': 79.8658, 'Al2O3': 101.9613,
            'Fe2O3': 159.6882, 'Cr2O3': 151.9904, 'FeO': 71.8444,
            'MnO': 70.9374, 'MgO': 40.3044, 'NiO': 74.6928,
            'CoO': 74.9326, 'CaO': 56.0774, 'Na2O': 61.9789,
            'K2O': 94.1960, 'P2O5': 141.9445, 'H2O': 18.0153,
            'CO2': 44.0095, 'SO3': 80.0632, 'Cl2O-1': 70.9045,
            'F2O-1': 37.9968
        }
        
        # 양이온 수
        self.cation_numbers = {
            'SiO2': 1, 'TiO2': 1, 'Al2O3': 2, 'Fe2O3': 2, 'Cr2O3': 2,
            'FeO': 1, 'MnO': 1, 'MgO': 1, 'NiO': 1, 'CoO': 1,
            'CaO': 1, 'Na2O': 2, 'K2O': 2, 'P2O5': 2, 'H2O': 2,
            'CO2': 1, 'SO3': 1, 'Cl2O-1': 2, 'F2O-1': 2
        }
        
        # 산소 수
        self.oxygen_numbers = {
            'SiO2': 2, 'TiO2': 2, 'Al2O3': 3, 'Fe2O3': 3, 'Cr2O3': 3,
            'FeO': 1, 'MnO': 1, 'MgO': 1, 'NiO': 1, 'CoO': 1,
            'CaO': 1, 'Na2O': 1, 'K2O': 1, 'P2O5': 5, 'H2O': 1,
            'CO2': 2, 'SO3': 3, 'Cl2O-1': -1, 'F2O-1': -1
        }
        
        # 광물 끝단 성분 데이터
        self.end_members = self.initialize_end_members()
        
    def initialize_end_members(self):
        """광물 끝단 성분 열역학 데이터 초기화"""
        # 간단한 예시 데이터 (실제로는 더 복잡한 데이터베이스 필요)
        return {
            'quartz': {
                'formula': 'SiO2',
                'H': -910700,  # J/mol
                'S': 41.46,    # J/mol/K
                'V': 22.688,   # cm3/mol
                'Cp_params': [46.94, 0.03435, -1129000, -241.8]  # Cp = a + bT + c/T^2 + d/T^0.5
            },
            'albite': {
                'formula': 'NaAlSi3O8',
                'H': -3935000,
                'S': 207.40,
                'V': 100.07,
                'Cp_params': [258.16, 0.05820, -2926200, -1226.6]
            },
            'anorthite': {
                'formula': 'CaAl2Si2O8',
                'H': -4228000,
                'S': 199.30,
                'V': 100.79,
                'Cp_params': [217.44, 0.10330, -3055300, -1165.6]
            },
            'forsterite': {
                'formula': 'Mg2SiO4',
                'H': -2174000,
                'S': 94.11,
                'V': 43.79,
                'Cp_params': [138.72, 0.02210, -1284300, -595.0]
            },
            'fayalite': {
                'formula': 'Fe2SiO4',
                'H': -1479000,
                'S': 150.79,
                'V': 46.39,
                'Cp_params': [138.72, 0.06480, -1258600, -595.0]
            },
            'enstatite': {
                'formula': 'MgSiO3',
                'H': -1545000,
                'S': 67.86,
                'V': 31.27,
                'Cp_params': [82.01, 0.01589, -714300, -308.2]
            },
            'ferrosilite': {
                'formula': 'FeSiO3',
                'H': -1195000,
                'S': 95.40,
                'V': 32.92,
                'Cp_params': [82.01, 0.03958, -677100, -308.2]
            }
        }
    
    def calculate_molar_composition(self, wt_composition: Dict[str, float]) -> Dict[str, float]:
        """중량% 조성을 몰 조성으로 변환"""
        molar_comp = {}
        total_moles = 0.0
        
        for oxide, wt_pct in wt_composition.items():
            if oxide in self.molecular_weights and wt_pct > 0:
                moles = wt_pct / self.molecular_weights[oxide]
                molar_comp[oxide] = moles
                total_moles += moles
        
        # 정규화
        if total_moles > 0:
            for oxide in molar_comp:
                molar_comp[oxide] /= total_moles
        
        return molar_comp
    
    def calculate_cation_fractions(self, molar_composition: Dict[str, float]) -> Dict[str, float]:
        """양이온 분율 계산"""
        cation_moles = {}
        total_cations = 0.0
        
        # 주요 양이온
        cation_map = {
            'SiO2': 'Si', 'TiO2': 'Ti', 'Al2O3': 'Al', 'Fe2O3': 'Fe3+',
            'FeO': 'Fe2+', 'MgO': 'Mg', 'CaO': 'Ca', 'Na2O': 'Na',
            'K2O': 'K', 'MnO': 'Mn', 'P2O5': 'P', 'Cr2O3': 'Cr'
        }
        
        for oxide, mol_frac in molar_composition.items():
            if oxide in cation_map:
                cation = cation_map[oxide]
                n_cations = self.cation_numbers.get(oxide, 1)
                cation_moles[cation] = mol_frac * n_cations
                total_cations += mol_frac * n_cations
        
        # 정규화
        if total_cations > 0:
            for cation in cation_moles:
                cation_moles[cation] /= total_cations
        
        return cation_moles
    
    def calculate_liquid_properties(self, composition: Dict[str, float], T: float, P: float) -> Dict[str, float]:
        """액체 상의 열역학적 성질 계산"""
        # Lange & Carmichael (1990) 부분 몰 부피 모델의 간단한 버전
        
        # 부분 몰 부피 (cm3/mol) - 1673K, 1 bar 기준
        V_ref = {
            'SiO2': 26.86, 'TiO2': 28.32, 'Al2O3': 37.42, 'Fe2O3': 44.28,
            'FeO': 13.97, 'MgO': 12.02, 'CaO': 16.90, 'Na2O': 29.03,
            'K2O': 46.34, 'H2O': 22.90
        }
        
        # 열팽창 계수 (dV/dT)
        dVdT = {
            'SiO2': 0.0, 'TiO2': 0.00724, 'Al2O3': 0.00262, 'Fe2O3': 0.00909,
            'FeO': 0.00292, 'MgO': 0.00327, 'CaO': 0.00374, 'Na2O': 0.00741,
            'K2O': 0.01191, 'H2O': 0.00950
        }
        
        # 압축률 (dV/dP)
        dVdP = {
            'SiO2': -0.000189, 'TiO2': -0.000231, 'Al2O3': -0.000226,
            'Fe2O3': -0.000253, 'FeO': -0.000045, 'MgO': -0.000027,
            'CaO': -0.000034, 'Na2O': -0.000240, 'K2O': -0.000675,
            'H2O': -0.000320
        }
        
        T_ref = 1673.0  # K
        P_ref = 1.0     # bar
        
        # 몰 조성 계산
        mol_comp = self.calculate_molar_composition(composition)
        
        # 부피 계산
        V_total = 0.0
        mass_total = 0.0
        
        for oxide, mol_frac in mol_comp.items():
            if oxide in V_ref:
                V_partial = V_ref[oxide]
                if oxide in dVdT:
                    V_partial += dVdT[oxide] * (T - T_ref)
                if oxide in dVdP:
                    V_partial += dVdP[oxide] * (P - P_ref)
                
                V_total += mol_frac * V_partial
                mass_total += mol_frac * self.molecular_weights[oxide]
        
        # 밀도 계산
        density = mass_total / V_total if V_total > 0 else 2.5
        
        # 점성도 계산 (Giordano et al., 2008 모델의 간단한 버전)
        # log η = A + B/(T - C)
        # 여기서는 SiO2 함량에 기반한 간단한 모델 사용
        SiO2_wt = composition.get('SiO2', 0)
        A = -4.55
        B = 4500 + 60 * SiO2_wt
        C = 600 - 2 * SiO2_wt
        
        log_viscosity = A + B / (T - C)
        viscosity = 10 ** log_viscosity
        
        # 열용량 (간단한 선형 모델)
        Cp = 80 + 0.02 * T + 0.5 * SiO2_wt
        
        # 엔탈피와 엔트로피 (간단한 근사)
        H = -1400000 - 100 * T  # J/mol
        S = 200 + 0.1 * T       # J/mol/K
        
        # 깁스 자유 에너지
        G = H - T * S
        
        return {
            'density': density,
            'viscosity': viscosity,
            'volume': V_total,
            'Cp': Cp,
            'H': H / 1000,  # kJ/mol로 변환
            'S': S,
            'G': G / 1000,  # kJ/mol로 변환
            'mass': mass_total
        }
    
    def calculate_saturation_surface(self, composition: Dict[str, float], T: float, P: float,
                                   phase: str) -> Tuple[bool, float]:
        """특정 상의 포화 여부와 친화도 계산"""
        # 간단한 포화 모델 (실제로는 복잡한 활동도 모델 필요)
        
        saturation_temps = {
            'quartz': {'base': 1000, 'SiO2_coeff': -5},
            'feldspar': {'base': 900, 'Al2O3_coeff': -10},
            'olivine': {'base': 1200, 'MgO_coeff': -20},
            'clinopyroxene': {'base': 1100, 'CaO_coeff': -15},
            'orthopyroxene': {'base': 1150, 'MgO_coeff': -18},
            'spinel': {'base': 1300, 'Al2O3_coeff': -8}
        }
        
        if phase not in saturation_temps:
            return False, 0.0
        
        # 포화 온도 계산
        sat_params = saturation_temps[phase]
        T_sat = sat_params['base']
        
        # 조성에 따른 보정
        for oxide, coeff_key in [('SiO2', 'SiO2_coeff'), ('Al2O3', 'Al2O3_coeff'),
                                 ('MgO', 'MgO_coeff'), ('CaO', 'CaO_coeff')]:
            if coeff_key in sat_params and oxide in composition:
                T_sat += sat_params[coeff_key] * composition[oxide]
        
        # 압력 보정 (간단한 선형 모델)
        T_sat += 0.03 * P  # 압력 증가시 포화 온도 상승
        
        # 과포화도 계산 (친화도의 근사)
        affinity = (T_sat - T) * 8.314 / 1000  # kJ/mol
        is_saturated = T <= T_sat
        
        return is_saturated, affinity
    
    def equilibrium_crystallization(self, bulk_composition: Dict[str, float], 
                                  T: float, P: float, fo2: float,
                                  allowed_phases: List[str]) -> Dict[str, Phase]:
        """평형 결정화 계산"""
        results = {}
        remaining_liquid = bulk_composition.copy()
        total_mass = sum(bulk_composition.values())
        
        # 액체 상 계산
        liquid_props = self.calculate_liquid_properties(remaining_liquid, T, P)
        results['liquid'] = Phase(
            name='liquid',
            mass=total_mass,
            composition=remaining_liquid,
            properties=liquid_props,
            formula=self.get_liquid_formula(remaining_liquid)
        )
        
        # 고체 상 검사
        crystallized_mass = 0.0
        
        for phase in allowed_phases:
            if phase == 'liquid':
                continue
                
            is_saturated, affinity = self.calculate_saturation_surface(remaining_liquid, T, P, phase)
            
            if is_saturated:
                # 간단한 결정화 모델
                crystal_fraction = min(0.1, abs(affinity) / 100)  # 최대 10%
                crystal_mass = total_mass * crystal_fraction
                
                # 광물 조성 계산 (간단한 예시)
                mineral_comp = self.calculate_mineral_composition(phase, remaining_liquid)
                
                results[phase] = Phase(
                    name=phase,
                    mass=crystal_mass,
                    composition=mineral_comp,
                    properties={'density': 3.0, 'affinity': affinity},
                    formula=self.get_mineral_formula(phase, mineral_comp)
                )
                
                crystallized_mass += crystal_mass
        
        # 액체 질량 업데이트
        results['liquid'].mass = total_mass - crystallized_mass
        
        # 전체 시스템 속성
        system_props = self.calculate_system_properties(results, T, P)
        results['system'] = Phase(
            name='system',
            mass=total_mass,
            composition=bulk_composition,
            properties=system_props
        )
        
        return results
    
    def fractional_crystallization(self, bulk_composition: Dict[str, float],
                                 T_start: float, T_end: float, T_step: float,
                                 P: float, fo2: float,
                                 allowed_phases: List[str]) -> List[Dict[str, Phase]]:
        """분별 결정화 계산"""
        results = []
        current_liquid = bulk_composition.copy()
        current_mass = sum(bulk_composition.values())
        
        temperatures = np.arange(T_start, T_end - T_step/2, -abs(T_step))
        
        for T in temperatures:
            # 현재 온도에서 평형 계산
            step_results = self.equilibrium_crystallization(
                current_liquid, T, P, fo2, allowed_phases
            )
            
            # 결정화된 광물 제거 (분별 결정화)
            total_crystal_mass = 0.0
            for phase_name, phase in step_results.items():
                if phase_name not in ['liquid', 'system']:
                    total_crystal_mass += phase.mass
                    
                    # 액체에서 원소 제거
                    for oxide in phase.composition:
                        if oxide in current_liquid:
                            removed = phase.mass * phase.composition[oxide] / 100
                            current_liquid[oxide] = max(0, current_liquid[oxide] - removed)
            
            # 액체 조성 재정규화
            liquid_sum = sum(current_liquid.values())
            if liquid_sum > 0:
                for oxide in current_liquid:
                    current_liquid[oxide] = current_liquid[oxide] / liquid_sum * 100
            
            current_mass -= total_crystal_mass
            results.append(step_results)
        
        return results
    
    def calculate_mineral_composition(self, phase: str, liquid_comp: Dict[str, float]) -> Dict[str, float]:
        """광물 조성 계산 (간단한 분배 계수 모델)"""
        # 실제로는 복잡한 열역학적 모델이 필요
        mineral_comp = {}
        
        if phase == 'olivine':
            # Mg-Fe 감람석
            Kd_Fe_Mg = 0.3  # Fe-Mg 분배계수
            mineral_comp['MgO'] = liquid_comp.get('MgO', 0) * 2.0
            mineral_comp['FeO'] = liquid_comp.get('FeO', 0) * Kd_Fe_Mg * 2.0
            mineral_comp['SiO2'] = 40.0
            
        elif phase == 'clinopyroxene':
            # 단사휘석
            mineral_comp['CaO'] = liquid_comp.get('CaO', 0) * 1.5
            mineral_comp['MgO'] = liquid_comp.get('MgO', 0) * 1.2
            mineral_comp['FeO'] = liquid_comp.get('FeO', 0) * 0.8
            mineral_comp['SiO2'] = 52.0
            mineral_comp['Al2O3'] = liquid_comp.get('Al2O3', 0) * 0.3
            
        elif phase == 'feldspar':
            # 장석
            mineral_comp['Al2O3'] = liquid_comp.get('Al2O3', 0) * 1.2
            mineral_comp['SiO2'] = 65.0
            mineral_comp['Na2O'] = liquid_comp.get('Na2O', 0) * 0.8
            mineral_comp['K2O'] = liquid_comp.get('K2O', 0) * 0.7
            mineral_comp['CaO'] = liquid_comp.get('CaO', 0) * 0.3
            
        elif phase == 'quartz':
            mineral_comp['SiO2'] = 100.0
            
        elif phase == 'orthopyroxene':
            # 사방휘석
            mineral_comp['MgO'] = liquid_comp.get('MgO', 0) * 1.5
            mineral_comp['FeO'] = liquid_comp.get('FeO', 0) * 1.0
            mineral_comp['SiO2'] = 55.0
            mineral_comp['Al2O3'] = liquid_comp.get('Al2O3', 0) * 0.2
        
        # 정규화
        total = sum(mineral_comp.values())
        if total > 0:
            for oxide in mineral_comp:
                mineral_comp[oxide] = mineral_comp[oxide] / total * 100
        
        return mineral_comp
    
    def get_mineral_formula(self, phase: str, composition: Dict[str, float]) -> str:
        """광물 화학식 계산"""
        cations = self.calculate_cation_fractions(self.calculate_molar_composition(composition))
        
        if phase == 'olivine':
            Mg = cations.get('Mg', 0)
            Fe = cations.get('Fe2+', 0)
            total = Mg + Fe
            if total > 0:
                return f"(Mg{Mg/total:.2f}Fe{Fe/total:.2f})₂SiO₄"
                
        elif phase == 'clinopyroxene':
            Ca = cations.get('Ca', 0)
            Mg = cations.get('Mg', 0)
            Fe = cations.get('Fe2+', 0)
            Al = cations.get('Al', 0)
            return f"Ca{Ca:.2f}(Mg{Mg:.2f}Fe{Fe:.2f})Si{2-Al:.2f}Al{Al:.2f}O₆"
            
        elif phase == 'feldspar':
            Na = cations.get('Na', 0)
            K = cations.get('K', 0)
            Ca = cations.get('Ca', 0)
            Al = cations.get('Al', 0)
            Si = cations.get('Si', 0)
            return f"(Na{Na:.2f}K{K:.2f}Ca{Ca:.2f})Al{Al:.2f}Si{Si:.2f}O₈"
            
        elif phase == 'quartz':
            return "SiO₂"
            
        elif phase == 'orthopyroxene':
            Mg = cations.get('Mg', 0)
            Fe = cations.get('Fe2+', 0)
            return f"(Mg{Mg:.2f}Fe{Fe:.2f})SiO₃"
            
        return phase
    
    def get_liquid_formula(self, composition: Dict[str, float]) -> str:
        """액체 조성을 간단한 형식으로 표현"""
        major_oxides = []
        for oxide in ['SiO2', 'Al2O3', 'FeO', 'MgO', 'CaO', 'Na2O', 'K2O', 'H2O']:
            if oxide in composition and composition[oxide] > 0.1:
                major_oxides.append(f"{oxide} {composition[oxide]:.1f}%")
        return ", ".join(major_oxides[:5])  # 상위 5개만 표시
    
    def calculate_system_properties(self, phases: Dict[str, Phase], T: float, P: float) -> Dict[str, float]:
        """전체 시스템의 열역학적 성질 계산"""
        total_mass = 0.0
        total_volume = 0.0
        total_H = 0.0
        total_S = 0.0
        total_Cp = 0.0
        
        for phase_name, phase in phases.items():
            if phase_name == 'system':
                continue
                
            total_mass += phase.mass
            
            if 'volume' in phase.properties:
                total_volume += phase.properties['volume'] * phase.mass / 100
            if 'H' in phase.properties:
                total_H += phase.properties['H'] * phase.mass / 100
            if 'S' in phase.properties:
                total_S += phase.properties['S'] * phase.mass / 100
            if 'Cp' in phase.properties:
                total_Cp += phase.properties['Cp'] * phase.mass / 100
        
        density = total_mass / total_volume if total_volume > 0 else 2.7
        G = total_H - T * total_S / 1000  # kJ
        
        return {
            'mass': total_mass,
            'volume': total_volume,
            'density': density,
            'H': total_H,
            'S': total_S,
            'G': G,
            'Cp': total_Cp,
            'T': T,
            'P': P
        }
    
    def calculate_oxygen_fugacity(self, T: float, P: float, buffer: str = "QFM") -> float:
        """산소 분압 계산"""
        # 버퍼별 log fO2 계산 (Frost, 1991)
        # log fO2 = A/T + B + C*(P-1)/T
        
        buffers = {
            'QFM': {'A': -24930, 'B': 8.29, 'C': 0.092},     # Quartz-Fayalite-Magnetite
            'NNO': {'A': -24930, 'B': 9.36, 'C': 0.092},     # Nickel-Nickel Oxide
            'IW': {'A': -27215, 'B': 6.57, 'C': 0.055},      # Iron-Wüstite
            'MW': {'A': -32807, 'B': 7.01, 'C': 0.035},      # Magnetite-Wüstite
            'FMQ': {'A': -25096, 'B': 8.74, 'C': 0.110}      # Fayalite-Magnetite-Quartz
        }
        
        if buffer not in buffers:
            return 0.0
        
        params = buffers[buffer]
        log_fo2 = params['A'] / T + params['B'] + params['C'] * (P - 1) / T
        
        return log_fo2
    
    def calculate_viscosity(self, composition: Dict[str, float], T: float, H2O_content: float = None) -> float:
        """멜트 점성도 계산 (Giordano et al., 2008)"""
        # VFT 방정식: log η = A + B/(T - C)
        # 여기서는 간단한 모델 사용
        
        if H2O_content is None:
            H2O_content = composition.get('H2O', 0)
        
        # 조성 파라미터
        SiO2 = composition.get('SiO2', 0)
        Al2O3 = composition.get('Al2O3', 0)
        alkalis = composition.get('Na2O', 0) + composition.get('K2O', 0)
        
        # VFT 파라미터 (조성 의존적)
        A = -4.55 + 0.034 * (SiO2 - 50)
        B = 4800 + 80 * (SiO2 - 50) - 600 * H2O_content
        C = 540 + 2.8 * (SiO2 - 50) + 14 * Al2O3 - 16 * alkalis
        
        log_viscosity = A + B / (T - C)
        
        return 10 ** log_viscosity
    
    def calculate_density(self, composition: Dict[str, float], T: float, P: float) -> float:
        """멜트 밀도 계산 (Lange & Carmichael, 1990)"""
        # 부분 몰 부피 방법
        
        # 1673K, 1 bar에서의 부분 몰 부피 (cm³/mol)
        V_ref = {
            'SiO2': 26.86, 'TiO2': 28.32, 'Al2O3': 37.42,
            'Fe2O3': 44.28, 'FeO': 13.97, 'MgO': 12.02,
            'CaO': 16.90, 'Na2O': 29.03, 'K2O': 46.34,
            'H2O': 22.90
        }
        
        # 열팽창 계수 (∂V/∂T)
        dVdT = {
            'SiO2': 0.0, 'TiO2': 0.00724, 'Al2O3': 0.00262,
            'Fe2O3': 0.00909, 'FeO': 0.00292, 'MgO': 0.00327,
            'CaO': 0.00374, 'Na2O': 0.00741, 'K2O': 0.01191,
            'H2O': 0.00950
        }
        
        # 압축률 (∂V/∂P) × 10⁴
        dVdP = {
            'SiO2': -1.89, 'TiO2': -2.31, 'Al2O3': -2.26,
            'Fe2O3': -2.53, 'FeO': -0.45, 'MgO': -0.27,
            'CaO': -0.34, 'Na2O': -2.40, 'K2O': -6.75,
            'H2O': -3.20
        }
        
        T_ref = 1673.0  # K
        P_ref = 1.0     # bar
        
        # 몰 조성으로 변환
        mol_comp = self.calculate_molar_composition(composition)
        
        # 총 부피와 질량 계산
        total_volume = 0.0
        total_mass = 0.0
        
        for oxide, mol_frac in mol_comp.items():
            if oxide in V_ref:
                # 부분 몰 부피 계산
                V = V_ref[oxide]
                
                if oxide in dVdT:
                    V += dVdT[oxide] * (T - T_ref)
                    
                if oxide in dVdP:
                    V += dVdP[oxide] * (P - P_ref) / 10000
                
                total_volume += mol_frac * V
                total_mass += mol_frac * self.molecular_weights[oxide]
        
        # 밀도 = 질량 / 부피
        density = total_mass / total_volume if total_volume > 0 else 2.5
        
        return density