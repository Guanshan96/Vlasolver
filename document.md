# Vlasolver使用指南

Vlasolver对无碰撞一维等离子体中的静电扰动进行动理学仿真。本文档将对Vlasolver使用的物理模型，配置文件结构及求解器调用方法进行说明。

## 物理模型

Vlasolver依据一维情形下的Vlasov-Poisson方程
$$\frac{\partial f}{\partial t}+v\frac{\partial f}{\partial x}+\frac{eE}{m}\frac{\partial f}{\partial v}=0$$
以及
$$\frac{\partial^2 \phi}{\partial x^2}=-\frac{\rho}{\epsilon_0}$$

