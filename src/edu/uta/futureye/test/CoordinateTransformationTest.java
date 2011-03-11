package edu.uta.futureye.test;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.GridLayout;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.BorderFactory;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.UIManager;

import edu.uta.futureye.function.Variable;

public class CoordinateTransformationTest {

	public Component createComponents() {
		final JLabel label = new JLabel("局部坐标 (r,s) => 物理坐标 (x,y)");
		final PanelDraw pane = new PanelDraw();
		pane.setBorder(BorderFactory.createEmptyBorder(50, // top
				50, // left
				700, // bottom
				700) // right
				);
		pane.setLayout(new GridLayout(0, 1)); // 单列多行
		pane.add(label);
		pane.addMouseListener(new MouseListener() {

			@Override
			public void mouseClicked(MouseEvent e) {
				// TODO Auto-generated method stub
				int x = e.getX();
				int y = e.getY();
				
				Variable p1 = new Variable();
				Variable p2 = new Variable();
				Variable p3 = new Variable();
				Variable p4 = new Variable();
				
				p1.set("r", (x-700)/100.0);
				p1.set("s", (y-200)/100.0);
				p2.set("r", (x-700)/100.0+0.25);
				p2.set("s", (y-200)/100.0);
				p3.set("r", (x-700)/100.0+0.25);
				p3.set("s", (y-200)/100.0+0.25);
				p4.set("r", (x-700)/100.0);
				p4.set("s", (y-200)/100.0+0.25);
				
				int x1 = (int)pane.fx.value(p1);
				int y1 = (int)pane.fy.value(p1);
				int x2 = (int)pane.fx.value(p2);
				int y2 = (int)pane.fy.value(p2);
				int x3 = (int)pane.fx.value(p3);
				int y3 = (int)pane.fy.value(p3);
				int x4 = (int)pane.fx.value(p4);
				int y4 = (int)pane.fy.value(p4);
				
				label.setText((x-700)/100.0+","+(y-200)/100.0+"  => "+"("+x1+","+y1+")");
				
				pane.getGraphics().drawRect(x,y,25,25);
				int [] xPoints = {x1,x2,x3,x4};
				int [] yPoints = {y1,y2,y3,y4};
				pane.getGraphics().drawPolygon(xPoints, yPoints, 4);
				


			}

			@Override
			public void mouseEntered(MouseEvent e) {
				// TODO Auto-generated method stub
				
			}

			@Override
			public void mouseExited(MouseEvent e) {
				// TODO Auto-generated method stub
				
			}

			@Override
			public void mousePressed(MouseEvent e) {
				// TODO Auto-generated method stub
				
			}

			@Override
			public void mouseReleased(MouseEvent e) {
				// TODO Auto-generated method stub
				
			}

		});		

        
		return pane;
	}

	public static void main(String[] args) {
		try {
			UIManager.setLookAndFeel(UIManager
					.getCrossPlatformLookAndFeelClassName());
			// 设置窗口风格
		} catch (Exception e) {
		}

		// 创建顶层容器并添加内容.
		JFrame frame = new JFrame("CoordinateTransformation Test");
		CoordinateTransformationTest app = new CoordinateTransformationTest();
		Component contents = app.createComponents();
		frame.getContentPane().add(contents, BorderLayout.CENTER);

		// 窗口设置结束，开始显示
		frame.addWindowListener(new WindowAdapter() {
			// 匿名类用于注册监听器
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});
		frame.pack();
		frame.setVisible(true);
	}

}
