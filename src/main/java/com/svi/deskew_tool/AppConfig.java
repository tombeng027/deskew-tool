package com.svi.deskew_tool;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public enum AppConfig {
	OUTPUT_PATH("OUTPUT_PATH"),
	IMAGE_INPUT_PATH("IMAGE_INPUT_PATH");
	
	private String value = "";
	private static Properties prop;

	private AppConfig(String value) {
		this.value = value;
	}

	public String value() {
		return prop.getProperty(value).trim();
	}

	public static void setContext(InputStream inputStream) {
		synchronized (inputStream) {
			if (prop == null) {
				try {
					prop = new Properties();
					prop.load(inputStream);
					inputStream.close();
				} catch (IOException e) {
					e.printStackTrace();
					try {
						inputStream.close();
					} catch (IOException e1) {
						e1.printStackTrace();
					}
				}
			}
		}
	}
}
