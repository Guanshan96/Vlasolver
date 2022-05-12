# Vlasolver使用指南

Vlasolver对无碰撞一维等离子体中的静电扰动进行动理学仿真。本文档将对Vlasolver使用的物理模型，配置文件结构及求解器调用方法进行说明。

## 物理模型

一维等离子体中的静电扰动由如下Vlasov-Poisson方程描述：
$$\frac{\partial f}{\partial t}+v\frac{\partial f}{\partial x}+\frac{eE}{m}\frac{\partial f}{\partial v}=0$$
以及
$$\frac{\partial^2 \phi}{\partial x^2}=-\frac{\rho}{\epsilon_0}.$$

Vlasov方程刻画了相空间体积元中代表点的密度沿方程特征线（即粒子在自洽静电场$E$作用下的运动轨迹）运动时密度保持不变这一图像。据此，Vlasolver利用半Lagrange法求解上述系统，即通过跟踪方程特征线确定分布函数$f$以及电荷密度$\rho$的演化。

## 配置文件结构

配置文件分为两部分，其一为求解器配置文件，其一为初始条件配置文件。

><Plasma>
  <Specie name="electron" save="true">
    <charge>-1</charge>
    <cmratio>-1</cmratio>
    <density>1</density>
    <xdistr number="1">
      <distr>@(x)1+0.01*cos(0.3*x)+0.01*cos(0.2*x)+0.01*cos(0.4*x)</distr>
    </xdistr>
    <vdistr number="1">
      <distr>@(vx)30*exp(-vx.^2/2)/(31*sqrt(2*pi))</distr>
    </vdistr>
  </Specie>
  <Specie name="electron" save="true">
    <charge>-1</charge>
    <cmratio>-1</cmratio>
    <density>1</density>
    <xdistr number="1">
      <distr>@(x)1</distr>
    </xdistr>
    <vdistr number="1">
      <distr>@(vx)exp(-(vx-6).^2/1)/(31*sqrt(pi))</distr>
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
