/*global Backbone, jQuery, _, ENTER_KEY */
var app = app || {};

$(function () {
	app.data = {};
	app.view = {}
	app.id = 0;
	app.data[app.id.toString()] = new app.MesmerData();
	app.curData = app.data[app.id.toString()];
	app.view[app.id.toString()] = new app.OverView({model: app.curData});
	app.id++;
});


