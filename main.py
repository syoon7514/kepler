import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import time

st.set_page_config(layout="wide")
st.title("\U0001F30C 태양계 행성의 케플러 법칙 시뮬레이터")

# 태양계 행성 데이터
planet_data = {
    "Mercury": {"a": 0.387, "e": 0.206, "T": 0.241},
    "Venus": {"a": 0.723, "e": 0.007, "T": 0.615},
    "Earth": {"a": 1.000, "e": 0.017, "T": 1.000},
    "Mars": {"a": 1.524, "e": 0.093, "T": 1.881},
    "Jupiter": {"a": 5.203, "e": 0.049, "T": 11.862},
    "Saturn": {"a": 9.537, "e": 0.056, "T": 29.457},
    "Uranus": {"a": 19.191, "e": 0.047, "T": 84.011},
    "Neptune": {"a": 30.07, "e": 0.009, "T": 164.8}
}

e_scale = 5  # 이심률 과장 배율
total_duration = 30  # 전체 애니메이션 시간 (초)
simulation_speed = 0.03  # 1프레임당 시간 (초)
total_frames = int(total_duration / simulation_speed)

# 행성 선택 UI
st.subheader("\U0001FA90 행성을 선택하세요")
cols = st.columns(len(planet_data))
selected_planet = None
for i, (name, _) in enumerate(planet_data.items()):
    if cols[i].button(name):
        selected_planet = name

if selected_planet:
    a = planet_data[selected_planet]["a"]
    e_real = planet_data[selected_planet]["e"]
    e = min(e_real * e_scale, 0.9)
    T = planet_data[selected_planet]["T"]

    # 공통 프레임 수를 가지면서 주기에 비례한 시간 간격 설정
    dt = T / total_frames

    GMsun = 4 * np.pi**2

    theta_all = np.linspace(0, 2*np.pi, 500)
    r_all = a * (1 - e**2) / (1 + e * np.cos(theta_all))
    x_orbit = r_all * np.cos(theta_all)
    y_orbit = r_all * np.sin(theta_all)

    positions, velocities, times, thetas, rs = [], [], [], [], []

    for step in range(total_frames):
        t = step * dt
        theta = 2 * np.pi * (t / T)
        r = a * (1 - e**2) / (1 + e * np.cos(theta))
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        v = np.sqrt(GMsun * (2/r - 1/a))

        positions.append((x, y))
        velocities.append(v * 30)  # 시각화를 위한 배율
        times.append(t)
        thetas.append(theta)
        rs.append(r)

    def sector_area(r1, r2, dtheta):
        return 0.5 * r1 * r2 * abs(dtheta)

    steps_20 = int(total_frames * 0.2)
    start_area = sum(
        sector_area(rs[i], rs[i+1], thetas[i+1] - thetas[i]) for i in range(steps_20-1)
    )
    end_area = sum(
        sector_area(rs[-i-2], rs[-i-1], thetas[-i-1] - thetas[-i-2]) for i in range(steps_20-1)
    )

    st.markdown(f"""
    **선택한 행성**: {selected_planet}  
    실제 이심률: {e_real:.3f} → 과장된 이심률: **{e:.3f}**  
    공전 반지름 a = {a:.3f} AU, 공전 주기 T = {T:.3f} 년

    ### \U0001F4D0 케플러 제2법칙: 면적 비교
    - 공전 초반 20% 부채꼴 면적: {start_area:.5f} AU²  
    - 공전 마지막 20% 부채꼴 면적: {end_area:.5f} AU²  
    👉 두 면적이 거의 동일함을 통해 **같은 시간 동안 같은 면적을 휩쓴다**는 법칙을 확인할 수 있어요.
    """)

    plot_area, graph_area = st.columns(2)

    for i in range(total_frames):
        x, y = positions[i]
        vx = -np.sin(thetas[i]) * velocities[i]
        vy = np.cos(thetas[i]) * velocities[i]

        # 공전 궤도 그래프
        fig1, ax1 = plt.subplots(figsize=(3.5, 3.5))
        ax1.plot(x_orbit, y_orbit, 'gray', lw=1)
        ax1.plot(0, 0, 'yo')
        ax1.plot(x, y, 'bo')
        ax1.quiver(x, y, vx, vy, color='red', scale=15, width=0.007)
        ax1.set_aspect('equal')
        ax1.set_xlim(-2*a, 2*a)
        ax1.set_ylim(-1.5*a, 1.5*a)
        ax1.set_xlabel("x (AU)")
        ax1.set_ylabel("y (AU)")
        ax1.set_title(f"{selected_planet} - Time: {times[i]:.2f} yr")
        ax1.grid(True)

        # 속도 그래프 (실시간 업데이트)
        fig2, ax2 = plt.subplots(figsize=(3.5, 3.5))
        ax2.plot(times[:i+1], velocities[:i+1], color='green')
        ax2.set_xlabel("Time (years)")
        ax2.set_ylabel("Orbital Speed (scaled km/s)")
        ax2.set_title("Orbital Speed vs Time")
        ax2.grid(True)

        with plot_area:
            st.pyplot(fig1)
        with graph_area:
            st.pyplot(fig2)

        time.sleep(simulation_speed)
else:
    st.info("행성을 선택하면 시뮬레이션이 시작됩니다.")
