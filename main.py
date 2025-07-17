import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import time

st.set_page_config(layout="wide")
st.title("🌞 태양계 케플러 법칙 시뮬레이터")

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

e_scale = 5

st.subheader("🌍 행성을 선택하세요")
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

    GMsun = 4 * np.pi**2
    dt = 0.03
    total_steps = int(T / dt)

    times = []
    thetas = []
    rs = []
    velocities = []

    for step in range(total_steps):
        t = step * dt
        theta = 2 * np.pi * (t / T)
        r = a * (1 - e**2) / (1 + e * np.cos(theta))
        v = np.sqrt(GMsun * (2/r - 1/a))

        times.append(t)
        thetas.append(theta)
        rs.append(r)
        velocities.append(v * 30)

    def sector_area(r1, r2, dtheta):
        return 0.5 * r1 * r2 * abs(dtheta)

    steps_20 = int(total_steps * 0.2)
    start_area_sector = sum(
        sector_area(rs[i], rs[i+1], thetas[i+1] - thetas[i]) for i in range(steps_20-1)
    )
    end_area_sector = sum(
        sector_area(rs[-i-2], rs[-i-1], thetas[-i-1] - thetas[-i-2]) for i in range(steps_20-1)
    )

    # 텍스트 먼저 출력
    st.markdown(f"""
    **선택한 행성**: {selected_planet}  
    실제 이심률: {e_real:.3f} → 과장된 이심률: **{e:.3f}**  
    공전 반지름 a = {a:.3f} AU, 공전 주기 T = {T:.3f} 년  

    ### 📐 케플러 제2법칙: 부채꼴 면적 계산  
    - **공전 초반 20% 부채꼴 면적**: {start_area_sector:.5f} AU²  
    - **공전 마지막 20% 부채꼴 면적**: {end_area_sector:.5f} AU²  
    👉 두 면적이 유사함을 통해 **면적 속도 일정성(케플러 제2법칙)**을 확인할 수 있습니다.
    """)

    # 궤도 전체 경로
    theta_all = np.linspace(0, 2*np.pi, 500)
    r_all = a * (1 - e**2) / (1 + e * np.cos(theta_all))
    x_orbit = r_all * np.cos(theta_all)
    y_orbit = r_all * np.sin(theta_all)

    plot_area = st.empty()
    graph_area = st.empty()

    for step in range(total_steps):
        theta = thetas[step]
        r = rs[step]
        t = times[step]
        v_scaled = velocities[step]

        x = r * np.cos(theta)
        y = r * np.sin(theta)
        vx = -v_scaled / 30 * np.sin(theta)
        vy = v_scaled / 30 * np.cos(theta)

        fig1, ax1 = plt.subplots(figsize=(6, 6))
        ax1.plot(x_orbit, y_orbit, 'gray', lw=1, label='Orbit Path')
        ax1.plot(0, 0, 'yo', label='Sun')
        ax1.plot(x, y, 'bo', label='Planet')
        ax1.quiver(x, y, vx, vy, color='red', scale=15, width=0.007, label='Velocity Vector')
        ax1.set_aspect('equal')
        ax1.set_xlim(-2*a, 2*a)
        ax1.set_ylim(-1.5*a, 1.5*a)
        ax1.set_xlabel("x (AU)")
        ax1.set_ylabel("y (AU)")
        ax1.set_title(f"{selected_planet} – Time = {t:.2f} yr")
        ax1.legend()
        ax1.grid(True)

        fig2, ax2 = plt.subplots()
        ax2.plot(times[:step+1], velocities[:step+1], color='green')
        ax2.set_xlabel("Time (years)")
        ax2.set_ylabel("Orbital Speed (scaled km/s)")
        ax2.set_title("Orbital Speed - Time")
        ax2.grid(True)

        with plot_area:
            st.pyplot(fig1)
        with graph_area:
            st.pyplot(fig2)

        time.sleep(0.03)
else:
    st.info("행성을 선택하면 시뮬레이션이 시작됩니다.")
