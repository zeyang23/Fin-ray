# Modeling of Fin-ray

This repo contains the codes of my undergraduate thesis (2021-06-01).

## Abstract

This article focuses on a kind of adaptive soft gripper—Fin-ray. Although Fin-ray has a wide range of applications in industrial production, its application scenarios are limited to open-loop systems. It is not integrated with soft sensors to form a closed-loop system, and is only used for adaptive grasping of objects, which cannot fully realize its potential as a soft manipulator. Due to the complex contact between Fin-ray and the object in actual working conditions and the large deflection, it is difficult to establish an accurate mechanical model to describe its "force-deflection" behaviors. The difficulty of modeling becomes the main obstacle to expanding its functions. In this context, this article proposes a large-deflection kinetostatic model of Fin-ray based on the principal axes decomposition of stiffness matrices, reflecting its "force-deflection" behaviors accurately. Furthermore, this article proposes the Fin-ray intrinsic sensing algorithm based on the information of the soft sensors on Fin-ray. The Fin-ray intrinsic sensing algorithm helps to establish a closed-loop feedback system, which enables Fin-ray to achieve dexterous manipulation and rich sensation at the same time.

## Some Results

### Roadmap of modeling

![image](https://user-images.githubusercontent.com/50078363/161173177-910dfaea-a8c8-4ca2-96da-ee8367c2fbd1.png)

### Validation of the theoretical model

![validation](https://user-images.githubusercontent.com/50078363/161171963-9bd3bd8e-60ca-40cf-bdc2-8e101d8970b2.gif)

### Intrinsic sensing and control of the grasping force

![抓取力大小的感知与控制](https://user-images.githubusercontent.com/50078363/161177093-b6e9797e-7ff8-433d-b9f7-373aba64b3f0.gif)

### Stiffness recognition

![物体刚度辨识](https://user-images.githubusercontent.com/50078363/161177139-66adf883-8928-4189-a28b-0f35089f0ff3.gif)

![image](https://user-images.githubusercontent.com/50078363/161177268-df881c1c-5069-4ca1-bc75-b31d9d0e9ead.png)

### Shape reconstruction

![物体形状辨识](https://user-images.githubusercontent.com/50078363/161176057-a9df9d72-12b0-417c-adf9-7e1af50fc887.gif)

### What's more

I used two methods to model the Fin-ray finger: Principal Axes Decomposition and Cosserat Rod Theory.

In the experiments, the model based on Principal Axes Decomposition is used.

However, the model based on Cosserat Rod Theory also has some advantages.

For example, it can deal with Fin-ray fingers with soft crossbeams, which is difficult for Principal Axes Decompoistion to solve.
It deals with the distributed load better than Principal Axes Decompoistion, which is an inherent discrete method.

![image](https://user-images.githubusercontent.com/50078363/161176799-b5008c52-3ccb-4b0c-a320-7bd894ad4293.png)

Besides, if we want to simulate the contact process more accurately, the line contact approximation is preferred rather than point contact.

Principal Axes Decomposition and Cosserat Rod Theory can be combined together to solve the line contact problem.

I wrote this in my undergraduate thesis, but I did not carry out experiments due to time constraints.

## Acknowledgement

I sincerely appreciate the help of my teammates: Zhijun Zhuang, Guangzhen Sun, and Erqi Tu.

A year has passed, I still miss the time we work together.

22-04-01
