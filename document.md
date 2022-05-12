# Vlasolver使用指南

Vlasolver对无碰撞一维等离子体中的静电扰动进行动理学仿真。本文档将对Vlasolver使用的物理模型，配置文件结构及求解器调用方法进行说明。

## 物理模型

一维等离子体中的静电扰动由如下Vlasov-Poisson方程描述：
$$\frac{\partial f}{\partial t}+v\frac{\partial f}{\partial x}+\frac{eE}{m}\frac{\partial f}{\partial v}=0$$
以及
$$\frac{\partial^2 \phi}{\partial x^2}=-\frac{\rho}{\epsilon_0}.$$

Vlasov方程刻画了相空间体积元中代表点的密度沿方程特征线（即粒子在自洽静电场$E$作用下的运动轨迹）运动时密度保持不变这一图像。据此，Vlasolver利用半Lagrange法求解上述系统，即通过跟踪方程特征线确定分布函数$f$以及电荷密度$\rho$的演化。不同于上述有量纲方程，模拟程序使用的是无量纲量，各物理参量的归一化按如下方式进行：

时间：
$$t=\omega_{pe}t'$$
长度:
$$x=\lambda_{De}x'$$
速度：
$$v=v_{Te}v'$$
电场强度：
$$E=\frac{m_e\omega_{pe}v_{Te}}{e}E'$$
电荷密度：
$$\rho=\frac{\epsilon_0m_e\omega_{pe}v_{Te}}{\lambda_{De}e}\rho'$$
其中$\omega_{pe}$为电子等离子体频率，$\lambda_{De}$为Debye长度，$v_{Te}$为电子热速度，$m_e$为电子静止质量。

## 求解器配置文件结构

以计算朗缪尔波静电衰变不稳定性为例，求解器配置文件形如
```
<Solver device="GPU"> <!--属性"device"决定代码在GPU或在CPU上运行-->
  <Grid>
    <xmax>4021.2385965949</xmax> <!--元素"xmax"决定系统位形空间的长度-->
    <xngrids>8192</xngrids>      <!--元素"xngrids"决定位形空间网格数-->
    <velGrid name="electron">    <!--元素"velGrid"决定速度空间网格参数，可按粒子成分配置-->
      <vmin>-12</vmin>           <!--元素"vmin"决定速度空间下界-->
      <vmax>20</vmax>            <!--元素"vmax"决定速度空间上界-->
      <vngrids>340</vngrids>     <!--元素"vngrids"决定速度空间网格数-->
    </velGrid>
    <velGrid name="proton">
      <vmin>-4</vmin>
      <vmax>4</vmax>
      <vngrids>256</vngrids>
    </velGrid>
  </Grid>
  <Temporal>
    <tstep>0.01</tstep>        <!--元素"tstep"决定时间步长-->
    <ntsteps>480000</ntsteps>  <!--元素"xmax"决定时间步数-->
  </Temporal>
  <Boundary>
    <type>PCHIP_periodic</type>       <!--元素"type"决定边界类型与插值方法-->
    <left name="Dirichlet">1</left>   <!--元素"left"决定非周期边界时Poisson方程的左边界条件-->
    <right name="Dirichlet">1</right> <!--元素"right"决定非周期边界时Poisson方程的右边界条件-->
  </Boundary>
  <Config>
    <filtration state="disable"> <!--速度空间条纹过滤。属性"state"决定启动与否-->
      <period>1000</period>      <!--元素"period"决定过滤时间周期-->
      <order>1</order>           <!--元素"order"决定滤波器阶数-->
      <width>5</width>           <!--元素"width"决定滤波器带宽-->
    </filtration>
    <external state="enable">         <!--外部电磁场。属性"state"决定启动与否-->
      <field type="global">           <!--元素"field"决定外场形式，属性"type"决定场的作用范围-->
        <magnetic>@(t, x)0</magnetic> <!--磁场（MATLAB函数句柄）-->
        <electric>@(t, x)0</electric> <!--电场（MATLAB函数句柄）-->
      </field>
      <field type="electron">
        <magnetic>@(t, x)0</magnetic>
        <electric>@(t, x)0.01*(1-exp(-t/30))*cos(0.15*x-t)</electric>
      </field>
      <field type="proton">
        <magnetic>@(t, x)0</magnetic>
        <electric>@(t, x)0</electric>
      </field>
    </external>
    <diagnostics state="enable"> <!--诊断参数。属性"state"决定启动与否-->
      <path>C:\Programmer\Vlasov-simulation_new\Diagnostics\Electrostatic1D_v1.m</path>
      <var>
        <name>Density</name> <!--诊断参数名称-->
        <rate>120</rate>     <!--诊断参数采样周期-->
      </var>
      <var>
        <name>Electric field</name>
        <rate>120</rate>
      </var>
      <var>
        <name>Distribution</name>
        <rate>12000</rate>
      </var>
    </diagnostics>
  </Config>
</Solver>
```
注释中标注了各个元素的意义与功能。

## 初始条件配置文件结构

初始条件配置文件形如
```
<Plasma>
  <Specie name="electron" save="true">     <!--粒子类型元素，属性"save"决定保存其分布函数与否-->
    <charge>-1</charge>                    <!--元素"charge"决定归一化电荷量（代数量）-->
    <cmratio>-1</cmratio>                  <!--元素"cmratio"决定归一化荷质比（代数量）-->
    <density>1</density>                   <!--元素"density"决定归一化密度-->
    <xdistr number="1">                    <!--位形空间分布函数，属性"number"决定分布函数数目（加和）-->
      <distr>@(x)1+0.01*cos(0.3*x)</distr> <!--位形空间分布函数（MATLAB函数句柄）-->
    </xdistr>
    <vdistr number="1">                            <!--速度空间分布函数，属性"number"决定分布函数数目（加和）---->
      <distr>@(vx)exp(-vx.^2/2)/sqrt(2*pi)</distr> <!--速度空间分布函数（MATLAB函数句柄）-->
    </vdistr>
  </Specie>
  <Specie name="proton" save="true">
    <charge>1</charge>
    <cmratio>0.1</cmratio>
    <density>1</density>
    <xdistr number="1">
      <distr>@(x)1+0.01*cos(0.3*x)+0.01*cos(0.2*x)+0.01*cos(0.4*x)</distr>
    </xdistr>
    <vdistr number="1">
      <distr>@(vx)exp(-vx.^2/2)/sqrt(2*pi)</distr>
    </vdistr>
  </Specie>
</Plasma>
```
